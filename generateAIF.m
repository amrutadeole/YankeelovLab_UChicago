%% Set-up in command window: 2821 voxels in ROI_1

imagesc(DCE_LD(:,:,1,15)); axis tight equal; %display image
ROI_1 = images.roi.Circle(gca,'Center',[272,151],'Radius',10); %draw ROI_1
bw = createMask(ROI_1); %create mask
[I,J] = find(bw); %find 1s: [row, column]
indices = [I J]; %store indices

%% ROI_1 AIF generate
% indices: contains [row, column] of slice 1 voxels in ROI_1
% DCE_LD_slice_no1: 400x400x30, slice 1 over time course

%% Iterate through indices and see if corresponding voxels meet AIF criteria
indicestf = ones(length(indices),1); %initialize 

for i = 1:length(indices)
    SIvoxel = DCE_LD_slice_no1(indices(i,1),indices(i,2),:); %will give you SI v. time for each voxel in ROI
    SI_max = max(SIvoxel);
    index_max = find(SIvoxel==SI_max); %find index of maximum SI value
   
    % If voxel does not meet criteria, mark as unusable
    if index_max > 5 %if maximum after 8th time point
        indicestf(i) = 0;
    end
    if SI_max < 10*std(SIvoxel(1:3))
        indicestf(i) = 0;
    end
    if mean(SIvoxel(20:30)) > 0.4*SI_max
        indicestf(i) = 0;
    end
   
end

%% Average SI of usable voxels to get AIF
goodvoxels = zeros(1,30); 
figure(3)
for i = 1:length(indicestf)
    if indicestf(i) == 1
        timecourse = DCE_LD_slice_no1(indices(i,1),indices(i,2),:);
        timecourse = reshape(timecourse,[1,30]);
        plot([1:1:30],timecourse)
        hold on
        ylim([0,5*10^6])
      
        %timecourse = timecourse';
        %goodvoxels = vertcat(goodvoxels, timecourse );
    end
end
goodvoxels(1,:) = []; %delete the first row of zeros

% Take average of all SI of usable voxels
AIF_ROI_1 = mean(goodvoxels);
        
%{
%_________________________________________________________________

%% Set-up: 2821 voxels in ROI_2
%{
imagesc(DCE_LD(:,:,1,1)); axis tight equal; %display image
ROI_2 = images.roi.Circle(gca,'Center',[277 161],'Radius',30); %draw ROI_2
bw = createMask(h); %create mask
[I,J] = find(bw); %find 1s: [row, column]
indices = [I J]; %store indices
%}

%% ROI_1 AIF generate
% indices: contains [row, column] of slice 1 voxels in ROI_1
% DCE_LD_slice_no1: 400x400x30, slice 1 over time course

%% Iterate through indices and see if corresponding voxels meet AIF criteria
indicesboolean = zeros(length(indices),1); %initialize 

for i = 1:length(indices)
    SIvoxel = DCE_LD_slice_no1(indices(i,1),indices(i,2),:); %will give you SI v. time for each voxel in ROI
    SI_max = max(SIvoxel);
    index_max = find(SIvoxel==SI_max); %find index of maximum SI value
    
    % If voxel does not meet criteria, mark as unusable
    if index_max > 8 %if maximum after 8th time point
        indicesboolean(i) = 0;
    end
    if SI_max < 20*std(SIvoxel(1:3))
        indicesboolean(i) = 0;
    end
    if mean(SIvoxel(20:30)) > 0.4*SI_max
        indicesboolean(i) = 0;
    end
end

%% Average SI of usable voxels to get AIF
goodvoxels = zeros(1,30); 
for i = 1:length(indicesboolean)
    if indicesboolean(i) == 1
        goodvoxels = vertcat(goodvoxels, DCE_LD_slice_no1(indices(i,1),indices(i,2),:));
    end
end
goodvoxels(1,:) = []; %delete the first row of zeros

% Take average of all SI of usable voxels
AIF_ROI_1 = mean(goodvoxels);
        
%}
        
        



        
        


