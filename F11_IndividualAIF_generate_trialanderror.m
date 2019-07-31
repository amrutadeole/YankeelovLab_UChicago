% 
% function F11_IndividualAIF_generate(subname)
% clear; 
% clc; 
% close all

% CASEname = 'AP049';
%subname = 'sub17';

%imagepath = '/Users/chengyuewu/Desktop/study&work/Lab/VesselTracking/Workspace/Data/orig/Cases/';
%matpath = '/Users/chengyuewu/Desktop/study&work/Lab/VesselTracking/Workspace/Data/Processed/';
%savepath = '/Users/chengyuewu/Desktop/study&work/Lab/VesselTracking/Workspace/Data/Processed/';

%load([matpath subname '_Artery_tracking_ITA.mat'])
%load([matpath subname '_DYNbr_ITA.mat'])
% load([origpath CASEname '.mat'],'DYN')


%arterynew = zeros(91,3);
%for i = 1:length(AIFSI_index)
%    arterynew(i,:) = (artery(AIFSI_index(i),:));
%end


FA = 30;
rel = 4.5;  %contrast agent relaxivity
TR_ht = 4.820/1000;

% R1 = 1./(T1_map/1000);
sina = sind(FA);
cosa = cosd(FA);
inj = 11;

T1_blood = 1.55;
T10 = T1_blood;
R10 = 1/T10;

reshaped_DCE_LD_array2 = permute(DCE_LD,[1,2,4,3]);
DCE_series = reshaped_DCE_LD_array2;
[ysize, xsize, numdyn, zsize] = size(DCE_series);
%%
SI = [];
indAIF_concentraiton = [];

if ~isempty(arterynew)
    for i = 1:size(arterynew,1)
        SI = [SI, permute(DCE_series(arterynew(i,1), arterynew(i,2), :, arterynew(i,3)), [3,1,2,4])];
    end

    SI = SI';

    %% signal -> average -> concentration
    S_pre = mean(mean(SI(:,1:inj),2));
%     ratio_0 = (1- exp(-TR_ht/T10))./(1- exp(-TR_ht/T10)*cosa); 
    S_0 = S_pre *(1-cosa*exp(-TR_ht/T10))/(sina*(1-exp(-TR_ht/T10)));
    % syms x
    % eqn = S_0*(sina*(1-exp(-TR_ht*(1/T10+rel*x)))) ./ (1-cosa*exp(-TR_ht*(1/T10+rel*x)))==SI_avg;

    % S_s = SI .* repmat(ratio_0, [size(artery,1), numDCE_series]) ./ repmat(S_pre, [1, numDCE_series]);
    % S = mean(S_s);

    SI_avg = mean(SI,1);
    S_norm = SI_avg/S_0;

%     indAIF_concentraiton = SI2C(S_norm,sina, cosa, TR_ht, T10,rel);
    indAIF_concentraiton = ((-log((S_norm-sina)./(S_norm*cosa-sina))/TR_ht - 1/T10)/rel);
    
    figure
    plot(1:length(indAIF_concentraiton), indAIF_concentraiton)
    xlabel('time points')
    ylabel('concentration (mM)')
    pause(0.1)
end
%% save data
%savename = [subname, '_individualAIF_ITA.mat'];
%save([savepath savename], 'indAIF_concentraiton','SI','T1_blood');


%%


%{
figure(3)
plot(3.4*(1:length(indAIF_concentraiton)), S_s')
xlabel('time(s)')
ylabel('SI / (S_0  sin\alpha)')

%
figure(2)
plot(3.4*(1:length(indAIF_concentraiton)), concentration')
xlabel('time(s)')
ylabel('concentration (mM)')

figure(4)
SS0 = S_s;
SS0(remove, :) = [];
plot(3.4*(1:length(indAIF_concentraiton)), SS0')
xlabel('time(s)')
ylabel('SI / (S_0  sin\alpha)')
%}

%%
% clear
% matpath = '/Users/chengyuewu/Desktop/study&work/Lab/VesselTracking/Workspace/Data/Processed/';
% subname = 'sub2';
% 
% load([matpath subname '_individualAIF.mat'])
% load([matpath subname '_DCE_seriesbr.mat'],'deltat');
% 
% initial = 0.1*[1,1,1,1,1];
% t = deltat * (1:(length(indAIF_concentraiton)));
% options = optimset('TolFun',1e-12,'Display','off','MaxIter',4000);
% 
% [para,~,resid] = lsqcurvefit('AIF_form',initial,t,indAIF_concentraiton, ...
%                               [],[],options);
%                           
% indAIF_fitted = AIF_form(para, t);
% 
% plot(indAIF_concentraiton)
% hold on
% plot(indAIF_fitted)
%         


