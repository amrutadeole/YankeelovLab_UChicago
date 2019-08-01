tic

% GENERATE DATA
%% Generate R(t) values to use as 'ydata'
T1map_sec = (T1map)./1000; %convert T10 map -> seconds
R1_sec = 1./(T1map_sec); %in 1/s
R1_sec(T1map_sec(:)==0) = 0;
s_t = double(reshaped_DCE_LD_array); %reshaped_DCE_LD_array: 320x320x24x30 instead of 400
TR = 4.820/1000; %repetition time in s
r1 = 4.5; %1/(mM*s)
flip = 10; %degrees
sina = sind(flip);
cosa = cosd(flip);

% Inherent Signal (S0 or M0), second 5 reps, before injection ~ index 12 (2
% min)
% Take mean of signal pre-injection(S) and solve for S0
s_inh = mean(s_t(:,:,:,6:10),4).*(1-exp(-TR*R1_sec).*cosd(flip))./((1-exp(-TR*R1_sec))*sind(flip)); %1:5 vandy data okay for S0
s_inh(isinf(s_inh)) = 0;
s_inh(isnan(s_inh)) = 0;

% R1 conversion
s_t_chosen = squeeze(s_t(160,160,12,:)); %reshaped 320x320x24x30 into 2457600x30 array

s_inh_chosen = s_inh(160,160,12); %reshaped 320x320x24 into 2457600x1 vector
for j = 1:30 %30 time points
    top = s_t_chosen(j)-s_inh_chosen*sina;
    bottom = s_t_chosen(j)*cosa-s_inh_chosen*sina;
    R_t_chosen(j) = (-1/TR)*log(top/bottom); %2457600x30 array 
end

R_t_chosen(isnan(R_t_chosen)) = 0;
R_t_chosen(isinf(R_t_chosen)) = 0;
R_t_chosen = abs(R_t_chosen);

%_______________________________________________________________________________
% FIT FOR KTRANS,VE
%% Initialize values
% Set no bounds
lb = [];
ub = [];

% Define AIF

timerange = [1:1:30];
timerange = TR.*timerange; %convert time to seconds

% Calculate initial guesses
ktrans_init = 0.1;
ve_init = 0.05;
%So=max(MFA_T1W_data(:));
params = [ktrans_init ve_init]; %params(1) = So*;  params(2) = ktrans; params(3) = ve

%% Fitting for So, ktrans, ve 
options = optimset('TolFun',1e-6,'Tolx',1e-6,'MaxIter',500,'Display','off');%,'algorithm','levenberg-marquardt');

param_out = zeros(1,2); %initialize param array for each T1 value: 2457600x3. column1=So, column2=ktrans, column3=ve
fits = zeros(1,30);
T10mapvoxel = T1map_sec(160,160,12); %get T10 value for chosen voxel
ydata = R_t_chosen;
[params, resnorm] = lsqcurvefit('t1_MFA_cf_KetyTofts',params,timerange,ydata,lb,ub,options,4.0668.*AIF_prostate_tenpercent,r1,T10mapvoxel);
calculated_ktrans = params(1);
calculated_ve = params(2);
param_out(1,:) = [calculated_ktrans calculated_ve];
tempfit = t1_MFA_cf_KetyTofts([params(1) params(2)],timerange,4.0668.*AIF_prostate_tenpercent,r1,T10mapvoxel);
fits(1,:) = tempfit;

figure(); plot(R_t_chosen,'.'); hold on; plot(fits);
%}
toc  