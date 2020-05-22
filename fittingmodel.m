clear all;
fdir1='/Users/Kalina/Desktop/Yankeelov Group/Data/ACRIN/DCEfiles/';
fdir2 ='/Users/Kalina/Desktop/Yankeelov Group/Data/ACRIN/Reg_DCE_images_vol/';
fdir3 = '/Users/Kalina/Desktop/Yankeelov Group/Data/ACRIN/muscle_R1/';
fdir_tumorR1 = '/Users/Kalina/Desktop/Yankeelov Group/Data/ACRIN/tumor_R1/';
fdir_new = '/Users/Kalina/Desktop/Yankeelov Group/Data/ACRIN/RR_ktrans_ve/test/';

sx = 256;
sy = sx;

% % original data
% ifreg = 0;
% fdir = 'H:\LX\Postdoc\BreastMR\ACRIN\DCEfiles\';
% sub_all=[15 20 22 27 84 88 143  439 453 460 464 495 508  714 720 724 732 733 867 882];

% % registered data
ifreg = 1;
fdir = '/Users/Kalina/Desktop/Yankeelov Group/Data/ACRIN/data_info.xlsx';
sub_all=[25, 58,141, 156,182,183,189,271,276,301,302,310,358, 516,529,530,632,651,716,718,723,725,770,793,815,824,832,839];

%sub_all=[20 25];
sub_all = [143];
[sz_all, nscan_all] = check_info(sub_all, fdir);
sn = 24;
global ktrans_RR ve_RR;

ktrans_RR = 0.15/60;
ve_RR = 0.12;
temp_res = 15;
TR = 0.1;
alpha = 90;
options=optimset('lsqcurvefit') ; 
options=optimset(options,'Display','off',...
    'TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',1e2*2,'MaxIter',2e3) ; 
%%
for ff = 1:length(sub_all)    
    sub = sub_all(ff)
    
    % % % load muscle R1 fit
%     load([fdir3 num2str(sub) '_R1_fit']);  
    load([fdir3 num2str(sub) '_R1_fit']);  
    
    if isnan(R1_fit) | isempty(R1_fit)
        continue;
    end 
    R1_RR = R1_fit;
    
    sz = sz_all(ff);
    nscan = nscan_all(ff);
    ktrans = zeros(sx,sy,sz);
    ve = zeros(sx,sy,sz);
    fits = zeros(sx*sy*sz,sn);
    t = [0:nscan-1]*temp_res;
        
    
    % % % load tumor R1 fit
%     load([fdir_tumorR1 num2str(sub) '_R1_fit']);  
    load([fdir_tumorR1 num2str(sub) '_R1_fit']);
    
    
    R1_TOI = R1_fit;
    pts = find(sum(R1_TOI,2)~=0); %tells us index of tumor voxels
    x0=[0.01, 0.1];
    xdata(:,1) = t;
    xdata(:,2) = R1_RR(:);

    for ii = 1:length(pts)
        jj = pts(ii);        
        ydata = R1_TOI(jj,:)';
        [x,resnorm] = lsqcurvefit(@func_RR_model,x0,xdata,ydata);
        ktrans(jj) = x(1);
        ve(jj) = x(2);
        tempfit = func_RR_model([x(1) x(2)], xdata, 1/R1_TOI(jj,1));
        fits(jj,:) = tempfit;
        
        %% 
%         x
%         R1_new = func_RR_model(x,xdata,1/R1_TOI(jj,1));
%         close all;
%         plot(ydata,'*');
%         hold on;
%         plot(R1_new, 'r.-');
%         pause;
        
    end
    clear xdata;
    
    ktrans=ktrans*60;
    pts = find(ktrans>5 | ktrans<=0.001 | ve<=0.001 | ve>=0.99);
%     ktrans(pts) = 0;
%     ve(pts) = 0;
    param_struct.ktrans = ktrans;
    param_struct.ve = ve;
    param_struct.patient = sprintf('%d',sub);
    param_struct.fits = fits;
%     save([fdir_new num2str(sub) '_para.mat '],'ktrans','ve','pts');
    save([fdir_new num2str(sub) '_para_test1.mat'],'param_struct');
end
