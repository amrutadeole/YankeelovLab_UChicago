for i = 1:size(arterynew2,1)
    plot(squeeze(DCE_LD(arterynew2(i,1),arterynew2(i,2),arterynew2(i,3),:)));
    hold on
end



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
fw = 1;

temp = zeros(1,30);
for i = 1:size(arterynew2,1)
    temp2 = squeeze(DCE_LD(arterynew2(i,1),arterynew2(i,2),arterynew2(i,3),:));
    temp2 = temp2';
    temp = temp + temp2;
end
avgS0 = temp/size(arterynew2,1); %outputs the average signal intensities of all the voxels in right gluteal muscle ROI

S_pre = mean(avgS0(1:11));
S_0 = S_pre *(1-cosa*exp(-TR_ht/T10))/(sina*(1-exp(-TR_ht/T10)));
S_norm = avgS0/S_0;
indAIF_concentration = ((-log((S_norm-sina)./(S_norm*cosa-sina))/TR_ht - 1/T10)/rel);
AIF_reference_gluteal = indAIF_concentration;
figure();plot(AIF_reference_gluteal) %outputs the AIF_gluteal_muscle
