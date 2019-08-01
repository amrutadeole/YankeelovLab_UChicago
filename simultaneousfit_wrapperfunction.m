

function [v_difference] = minimize_difference(guessvalues)
    AIF_prostate = [0.0150162856012648,0.0124456740013982,-0.0136446641280534,-0.00240130489801187,-0.00141583650913729,0.00352447644927720,0.0128797620562517,-0.00974867392689001,-0.0125073618573586,-0.00498463406867468,0.000988528235708028,0.00170670965564453,0.109991896948052,0.104698715021447,0.0440443675854020,0.0452625550214323,0.0498938911549105,0.0280469166480799,0.0267477018085782,0.0323492076365658,0.0294164260876530,0.0435361880571498,0.0347776505584428,0.0369105204427508,0.0288929814759390,0.0323142786744611,0.0374002204250499,0.0400685751375033,0.0225221544746800,0.0304024658787544];
    load('T1map_sec')
    truevalues = [0.075,0.35,0.5];
    r1 = 4.5; voxelindex = [160;160;12];
    timerange_FXL = [1:1:30];
    timerange_FXR = [1:1:65];
    %guessvalues = [0.06;0.2;0.4];
    Cp_FXL = 5.3084.*AIF_prostate;
    Cp_FXR = [0.0068,0.0019,0.0018,0.0012,0.0018,-0.0134,-0.0077,-0.0064,-0.0013,0.0052,0.0045,0.0845,0.1051,0.1848,0.1708,0.1806,0.1752,0.1693,0.1686,0.1626,0.1659,0.1668,0.1655,0.1608,0.1584,0.1561,0.1548,0.1593,0.1554,0.1523,0.1421,0.1374,0.1418,0.1374,0.141,0.14,0.137,0.1389,0.1341,0.1369,0.13,0.1318,0.1314,0.1311,0.1237,0.1262,0.1259,0.1279,0.1299,0.1286,0.1239,0.1285,0.129,0.1303,0.1361,0.1248,0.1234,0.1132,0.1253,0.1297,0.1297,0.1297,0.1297,0.1297,0.1297]; %NEEDS TO BE FINISHED
    TR = 4.820/1000;  
    fw = 1; 
    flip = 10;
    sina = sind(flip);
cosa = cosd(flip);
R1map_sec = 1./(T1map_sec); %in 1/s
    %R1 = 1./(T1map);
R1map_sec(T1map_sec(:)==0) = 0;
    %R1(T1map(:)==0) = 0;
s_t_LD = double(reshaped_DCE_LD_array); %reshaped_DCE_LD_array: 320x320x24x30 instead of 400
s_t_SD = double(reshaped_DCE_SD_array); %reshaped_DCE_SD_array: 320x320x24x30 instead of 400

% Inherent Signal (S0): take mean of signal pre-injection(S) and solve for S0
%s_inh = mean(s_t(:,:,:,6:10),4).*(1-exp(-TR*R1_sec).*cosd(flip))./((1-exp(-TR*R1_sec))*sind(flip)); %1:5 vandy data okay for S0
    s_inh_LD = mean(s_t_LD(:,:,:,6:10),4).*(1-exp(-TR*R1map_sec).*cosd(flip))./((1-exp(-TR*R1map_sec))*sind(flip)); %1:5 vandy data okay for S0
    s_inh_SD = mean(s_t_LD(:,:,:,6:10),4).*(1-exp(-TR*R1map_sec).*cosd(flip))./((1-exp(-TR*R1map_sec))*sind(flip)); %1:5 vandy data okay for S0

s_inh_LD(isinf(s_inh_LD)) = 0;
s_inh_LD(isnan(s_inh_LD)) = 0;
s_inh_SD(isinf(s_inh_LD)) = 0;
s_inh_SD(isnan(s_inh_LD)) = 0;
    %T1map_sec_reshaped = T1map_sec(:); %T10 values 
    %R1map_sec_reshaped = 1./T1map_sec_reshaped; %R10 values

    LD_SI_vector = LDtermfunction([truevalues(1),truevalues(2)],timerange_FXL,Cp_FXL,voxelindex,r1,flip,TR,s_inh_LD,R1map_sec);
    FXL_SI_vector = LDtermfunction([guessvalues(1),guessvalues(2)],timerange_FXL,Cp_FXL,voxelindex,r1,flip,TR,s_inh_LD,R1map_sec);
    a = LD_SI_vector - FXL_SI_vector;
    %x0st_LD = [0.1 0.2];
    %ydata_LD = LD_SI_vector;
    %[x,resnorm] = lsqcurvefit(@LDtermfunction,x0st_LD,timerange_FXL,ydata_LD,[],[],options,Cp_FXL,voxelindex,r1,flip,TR,s_inh_LD,R1map_sec);
    %LD_parameters = x
    %tempfit_LD = LDtmfunction([x(1) x(2)],timerange_FXL,Cp_FXL,voxelindex,r1,flip,TR,s_inh_LD,R1map_sec);
    %fits_LD = tempfit_LD;
    
    %SD_SI_vector = SDtermfunction([ktrans,ve,T_i],timerange_FXR,Cp_FXR,voxelindex,r1,flip,TR,fw,s_inh,R1map_sec);
    SD_SI_vector = SIequationwithnoise_FXR_7_30([truevalues(1),truevalues(2),truevalues(3)],timerange_FXR,Cp_FXR,r1,fw,flip,TR,s_inh_SD,R1map_sec,voxelindex);
    FXR_SI_vector = SIequationwithnoise_FXR_7_30([guessvalues(1),guessvalues(2),guessvalues(3)],timerange_FXR,Cp_FXR,r1,fw,flip,TR,s_inh_SD,R1map_sec,voxelindex);
    %figure();plot(SD_SI_vector);hold on;plot(FXR_SI_vector)
    b = SD_SI_vector - FXR_SI_vector;
    v_difference = [a,b];
    plot(timerange_FXR,FXR_SI_vector,'.k',timerange_FXR,SD_SI_vector,'-r',...
        timerange_FXL,FXL_SI_vector,'.b',timerange_FXL,LD_SI_vector,'-b');
    drawnow
  
end

    %x0st_SD = [0.1 0.2 0.4];time
    
    %ydata_SD = SD_SI_vector;
    %[x,resnorm] = lsqcurvefit(@SIequationwithnoise_FXR_7_30,x0st_SD,timerange_FXR,ydata_SD,[],[],options,Cp_FXR,r1,fw,flip,TR,s_inh_voxel,0,R1map_sec,voxelindex);
    %SD_parameters = x
    %tempfit_SD = SIequationwithnoise_FXR_7_30([x(1) x(2) x(3)],timerange_FXR,Cp_FXR,r1,fw,flip,TR,s_inh_voxel,0,R1map_sec,voxelindex);
    %fits_SD = tempfit_SD;
    
    %vector_FXL = fxltermfunction(timerange_FXL,Cp_FXL,[ktrans,ve],voxelindex,r1,flip,TR,s_inh_LD,R1map_sec);
    %vector_FXL = SIequationwithnoise_FXL_7_29([ktrans,ve],timerange_FXL,Cp_FXL,flip,TR,R1map_sec,r1,s_inh,voxelindex)
    %vector_FXL = vector_FXL';
    %vector_FXR = fxrtermfunction(timerange_FXR,Cp_FXR,[ktrans;ve;T_i],voxelindex,r1,flip,TR,fw,s_inh_SD,R1map_sec);
    %vector_FXR = vector_FXR';
    %vector_DCE_FXL = squeeze(reshaped_DCE_LD_array(voxelindex(1),voxelindex(2),voxelindex(3),:));
    %figure();plot(vector_FXL,'.');title('FXL generated');
    %figure();plot(vector_DCE_FXL);title('FXL DCE');
    %vector_DCE_FXR = squeeze(reshaped_DCE_SD_array(voxelindex(1),voxelindex(2),voxelindex(3),:));
    %figure();plot(vector_FXR,'.');title('FXR generated');
    %figure();plot(vector_DCE_FXR);title('FXR DCE');
    %v_LHS = [vector_FXL;vector_FXR];
    %v_RHS = [vector_DCE_FXL;vector_DCE_FXR];
    %v_difference = v_LHS - v_RHS;
%end
