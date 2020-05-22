% Forward Fitting Curve

tic
load('storerandom')
% Define assigned variables
assigned_ve = [0.05:0.05:0.5];
assigned_ktrans = [0.1:0.1:1.0];
assigned_ktrans = assigned_ktrans./60;

% Create random value array 1x66
%for i = 1:66
%    storerandom(i) = normrnd(0,1);
%end

% Initialize variables
S0 = 364.45; %chosen value from S0 map
Rel = 4.5; %mM^-1*s^-1
R10 = (1/1.5); %s^-1
flip = 5; %degrees
TR = 0.01; %s

rows = size(c_toi_values,1); 
columns = size(c_toi_values,2);

% Convert flip angle to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

% Initialize R1_vox matrix
R1_vox = zeros(rows,columns);

% Forward fitting 
parfor i = 1:length(assigned_ve) %i = "vertical slice" number
    x = [assigned_ktrans(i) assigned_ve(i)]; %a new vertical "slice" of the cube
    for k = 1:100 %k = # of runs per "vertical slice"
        SI_FXL = SIequationwithnoise_FXL(x,t,Cp,flip,TR,R10,Rel,S0,storerandom)
        SI_FXL = SI_FXL';
        SI_3d_FXL(:,k,i) = SI_FXL; %add vector to 3d array
    end
end


% LSQ fit
x0st=[0.1 0.2]; %initial guesses
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',150,'MaxFunEval',2000);
fits_FXL = zeros(66,100,10); %initialize fits 3d array
parfor ii = 1:10 %for each "vertical slice"
    for kk = 1:100 %for each run
        ydata = SI_3d_FXL(:,kk,ii);
        [x,resnorm] = lsqcurvefit('amrutaToftsFunSI',x0st,t_new,ydata,[],[],options,AIF_pop,flip,TR,R10,r1,S0);
        calculated_ktrans_FXL(kk,ii) = x(1);
        calculated_ve_FXL(kk,ii) = x(2);
        tempfit = amrutaToftsFunSI([x(1) x(2)],t_new,AIF_pop,flip,TR,R10,r1,S0);
        fits_FXL(:,kk,ii) = tempfit;
    end
end

toc

%{

% Calculate means of each parameter
for i = 1:10
    mean_ktrans_FXL(i) = sum(calculated_ktrans_FXL(:,i))./100;
    mean_ve_FXL(i) = sum(calculated_ve_FXL(:,i))./100;
    CI_ktrans_FXL_lower(i) = mean_ktrans_FXL(i) - 1.96*std(calculated_ktrans_FXL(:,i))/10;
    CI_ktrans_FXL_upper(i) = mean_ktrans_FXL(i) + 1.96*std(calculated_ktrans_FXL(:,i))/10;
    CI_ve_FXL_lower(i) = mean_ve_FXL(i) - 1.96*std(calculated_ve_FXL(:,i))/10;
    CI_ve_FXL_upper(i) = mean_ve_FXL(i) + 1.96*std(calculated_ve_FXL(:,i))/10;
end

fprintf('\t%1.4f\t\t   \t%1.2f\t\t\t%1.4f\t\t%1.2f %1.4f %1.4f \t\t %1.4f %f %2.3f\n',iter,eval,ea,unnormalized,normalized); %print calculated row on table
        
assigned_ktrans = assigned_ktrans';
assigned_ve = assigned_ve';
mean_ktrans_FXL = mean_ktrans_FXL';
mean_ve_FXL = mean_ve_FXL';
CI_ktrans_FXL_lower = CI_ktrans_FXL_lower';
CI_ktrans_FXL_upper = CI_ktrans_FXL_upper';
CI_ve_FXL_lower = CI_ve_FXL_lower';
CI_ve_FXL_upper = CI_ve_FXL_upper';



format short
%T = table(assigned_ktrans, assigned_ve, calculated_ktrans, calculated_ve)
Table_FXL = table(assigned_ktrans, assigned_ve, mean_ktrans_FXL, mean_ve_FXL, CI_ktrans_FXL_lower, CI_ktrans_FXL_upper, CI_ve_FXL_lower, CI_ve_FXL_upper)

toc 


%Calculate x(3) aka alpha value
x0st=[0.1 0.2 0.1]
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',150,'MaxFunEval',2000);
fits_FXL_alpha = zeros(66,100,10);
parfor ii = 1:10
    for kk = 1:100
        ydata = SI_3d_FXL(:,kk,ii);
        [x,resnorm] = lsqcurvefit('SIequationwithnoise_anda_FXL',x0st,t_new,ydata,[],[],options,AIF_pop,flip,TR,R10,r1,S0,storerandom);
        calculated_ktrans_FXL(kk,ii) = x(1);
        calculated_ve_FXL(kk,ii) = x(2);
        calculated_alpha_FXL(kk,ii) = x(3);
        tempfit = SIequationwithnoise_anda_FXL([x(1) x(2) x(3)],t_new,AIF_pop,flip,TR,R10,r1,S0,storerandom);
        fits_FXL_alpha(:,kk,ii) = tempfit;    
    end
end

toc

%calculate means of each parameter
mean_ktrans_FXL = sum(calculated_ktrans_FXL)./100;
mean_ve_FXL = sum(calculated_ve_FXL)./100;
mean_Ti_FXL = sum(calculated_Ti_FXL)./100;


figure(2)
plot(t_new,SI_3d_FXL(:,1,1),'.');
hold on
plot(t_new,fits_FXL(:,1,1));

toc 

% Plot data and fit on graphs
for jj = 1:10
    figure(jj)
    plot(t_new,S(:,jj), '.');
    xlabel('Time(sec)');
    ylabel('SI(arb. units');
    str = sprintf('SI Vs. Time Curve with ktrans = %4.4f and ve = %4.2f', assigned_ktrans(jj), assigned_ve(jj)); title(str);
    hold on
    plot(t_new,fits(jj,:));
    ylim([20,35]);
    hold off
end

% Plot family of curves on same graph
for jj = 1:10
    figure(11)
    plot(t_new,S(:,jj));
    xlabel('Time(sec)');
    ylabel('SI(arb. units)');
    title('Family of Curves');
    hold on
end

 
%}


