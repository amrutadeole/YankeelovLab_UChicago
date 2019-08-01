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


FA = 10;
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
%% Plot Top 10% AIFs
SI = [];
indAIF_concentraiton = [];

if ~isempty(chosen_artery10)
    for i = 1:size(chosen_artery10,1)
        SI = [SI, permute(DCE_series(chosen_artery10(i,1), chosen_artery10(i,2), :, chosen_artery10(i,3)), [3,1,2,4])];
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
    AIF_prostate_tenpercent = indAIF_concentraiton;
    figure (1)
    plot(1:length(indAIF_concentraiton), indAIF_concentraiton)
    xlabel('time points')
    ylabel('concentration (mM)')
    title('Top 10% AIF')
    pause(0.1)
end
%% Plot Top 15% AIFs
SI = [];
indAIF_concentraiton = [];

if ~isempty(chosen_artery15)
    for i = 1:size(chosen_artery15,1)
        SI = [SI, permute(DCE_series(chosen_artery15(i,1), chosen_artery15(i,2), :, chosen_artery15(i,3)), [3,1,2,4])];
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
    
    figure (2)
    plot(1:length(indAIF_concentraiton), indAIF_concentraiton)
    xlabel('time points')
    ylabel('concentration (mM)')
    title('Top 15% AIF')
    pause(0.1)
end
%% Plot Top 20% AIFs
SI = [];
indAIF_concentraiton = [];

if ~isempty(chosen_artery20)
    for i = 1:size(chosen_artery20,1)
        SI = [SI, permute(DCE_series(chosen_artery20(i,1), chosen_artery20(i,2), :, chosen_artery20(i,3)), [3,1,2,4])];
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
    
    figure (3)
    plot(1:length(indAIF_concentraiton), indAIF_concentraiton)
    xlabel('time points')
    ylabel('concentration (mM)')
    title('Top 20% AIF')
    pause(0.1)
end
%% Plot Top 25% AIFs
SI = [];
indAIF_concentraiton = [];

if ~isempty(chosen_artery25)
    for i = 1:size(chosen_artery25,1)
        SI = [SI, permute(DCE_series(chosen_artery25(i,1), chosen_artery25(i,2), :, chosen_artery25(i,3)), [3,1,2,4])];
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
    
    figure (4)
    plot(1:length(indAIF_concentraiton), indAIF_concentraiton)
    xlabel('time points')
    ylabel('concentration (mM)')
    title('Top 25% AIF')
    pause(0.1)
end
