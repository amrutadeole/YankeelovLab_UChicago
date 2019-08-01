

%% Initial values
enh0 = 12; %% for fast series  

%%% subtraction
DCE_LD_permuted = permute(DCE_LD,[1,2,4,3]); 
[ysize, xsize, numdyn, zsize] = size(DCE_LD_permuted); %DCE_LD: [y,x,time,# of slices]
injtime = 11; %injection time x/30

SUB = DCE_LD_permuted(:,:,(injtime+1):(end-injtime),:) - repmat(mean(DCE_LD_permuted(:,:,1:injtime,:), 3), [1 1 numdyn-2*injtime 1]); %400x400x8x24
pre_std = std(DCE_LD_permuted(:,:,1:injtime,:),[],3); %1x1

post_sub = permute(mean(SUB, 3), [1 2 4 3]); %400x400x24

%% Select seed point: left
P0_left = [166,120,12]; %set seed point [166,120], slice 12 of left artery
slice0 = P0_left(3); %slice of seed point

% Generate the initial kernel
[Y, X] = meshgrid((P0_left(1)-1):(P0_left(1)+1), (P0_left(2)-1):(P0_left(2)+1));
kernel0 = zeros(9,3); %9x3 matrix: columns are y, x, and slice#
kernel0(:,1) = Y(:);
kernel0(:,2) = X(:);
kernel0(:,3) = slice0;
artery_left = kernel0;

% tracking from one seed point
% Note: adjustment -- based on time-course correlation

% up  
slice = slice0;
P = P0_left;
kernel_intensity = post_sub((P0_left(1)-1):(P0_left(1)+1), (P0_left(2)-1):(P0_left(2)+1), slice);
kernel_intensity = kernel_intensity(:);
kernel_timecourse = mean(mean(DCE_LD_permuted((P0_left(1)-1):(P0_left(1)+1),(P0_left(2)-1):(P0_left(2)+1),:, slice)));
kernel_timecourse = kernel_timecourse(:);

for i = 1:(zsize-slice0) %for i = 1:# slices-slice of seed point
    slice = slice+1;
    [Y, X] = meshgrid((P(1)-2):(P(1)+2), (P(2)-2):(P(2)+2));
    ROI = zeros(25,3); %make 5x5 kernel for each slice
    ROI(:,1) = Y(:);
    ROI(:,2) = X(:);
    ROI(:,3) = slice;
    
    cc = 0;
    for j = 1:size(ROI,1) %for j = 1:25
        ptest = ROI(j,:);
        [y, x] = meshgrid((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1));
        window = zeros(9,2);
        window(:,1) = y(:);
        window(:,2) = x(:);
        window(:,3) = slice; 
        
        window_intensity = post_sub((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1), slice);
        window_intensity = window_intensity(:);
        window_timecourse = mean(mean(DCE_LD_permuted((ptest(1)-1):(ptest(1)+1),(ptest(2)-1):(ptest(2)+1),:, slice)));
        window_timecourse = window_timecourse(:);

        %calculate coorelation coefficient of SI between kernel and each
        %local window
        cc_update = corr(window_intensity, kernel_intensity);
        cc_update = corr(window_timecourse, kernel_timecourse);        
        if cc_update > cc
            P = ptest;
            kernel = window;
            cc = cc_update;
        end
    end
    
    kernel_intensity = post_sub((P(1)-1):(P(1)+1), (P(2)-1):(P(2)+1), slice);
    kernel_intensity = kernel_intensity(:);
    kernel_timecourse = mean(mean(DCE_LD_permuted((P(1)-1):(P(1)+1),(P(2)-1):(P(2)+1),:, slice)));
    kernel_timecourse = kernel_timecourse(:);
    
    artery_left = [artery_left; kernel];
end
    
%%%%% down
slice = slice0;
P = P0_left;
kernel_intensity = post_sub((P0_left(1)-1):(P0_left(1)+1), (P0_left(2)-1):(P0_left(2)+1), slice);
kernel_intensity = kernel_intensity(:);
kernel_timecourse = mean(mean(DCE_LD_permuted((P0_left(1)-1):(P0_left(1)+1),(P0_left(2)-1):(P0_left(2)+1),:, slice)));
kernel_timecourse = kernel_timecourse(:);

for i = 1:(slice0-1)
    slice = slice-1;
    [Y, X] = meshgrid((P(1)-2):(P(1)+2), (P(2)-2):(P(2)+2));
    ROI = zeros(25,3);
    ROI(:,1) = Y(:);
    ROI(:,2) = X(:);
    ROI(:,3) = slice;
    
    cc = 0;
    for j = 1:size(ROI,1)
        ptest = ROI(j,:);
        [y, x] = meshgrid((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1));
        window = zeros(9,2);
        window(:,1) = y(:);
        window(:,2) = x(:);
        window(:,3) = slice; 
        
        window_intensity = post_sub((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1), slice);
        window_intensity = window_intensity(:);
        window_timecourse = mean(mean(DCE_LD_permuted((ptest(1)-1):(ptest(1)+1),(ptest(2)-1):(ptest(2)+1),:, slice)));
        window_timecourse = window_timecourse(:);
        
        cc_update = corr(window_intensity, kernel_intensity);
        cc_update = corr(window_timecourse, kernel_timecourse);
        if cc_update > cc
            P = ptest;
            kernel = window;
            cc = cc_update;
        end
    end
    
    kernel_intensity = post_sub((P(1)-1):(P(1)+1), (P(2)-1):(P(2)+1), slice);
    kernel_intensity = kernel_intensity(:);
    kernel_timecourse = mean(mean(DCE_LD_permuted((P(1)-1):(P(1)+1),(P(2)-1):(P(2)+1),:, slice)));
    kernel_timecourse = kernel_timecourse(:);
%   
    artery_left = [artery_left; kernel];
end

%% Select seed point: right
P0_right = [166,277,12]; %set seed point [166,277], slice 12 of right artery
slice0 = P0_right(3);

%%%%%% generate the initial kernel
[Y, X] = meshgrid((P0_right(1)-1):(P0_right(1)+1), (P0_right(2)-1):(P0_right(2)+1));
kernel0 = zeros(9,3);
kernel0(:,1) = Y(:);
kernel0(:,2) = X(:);
kernel0(:,3) = slice0;
artery_right = kernel0;

%%% tracking from one seed point
%%%% Note: adjustment -- based on time-course correlation

%%%%% up  
slice = slice0;
P = P0_right;
kernel_intensity = post_sub((P0_right(1)-1):(P0_right(1)+1), (P0_right(2)-1):(P0_right(2)+1), slice);
kernel_intensity = kernel_intensity(:);
kernel_timecourse = mean(mean(DCE_LD_permuted((P0_right(1)-1):(P0_right(1)+1),(P0_right(2)-1):(P0_right(2)+1),:, slice)));
kernel_timecourse = kernel_timecourse(:);

for i = 1:(zsize-slice0)
    slice = slice+1;
    [Y, X] = meshgrid((P(1)-2):(P(1)+2), (P(2)-2):(P(2)+2));
    ROI = zeros(25,3);
    ROI(:,1) = Y(:);
    ROI(:,2) = X(:);
    ROI(:,3) = slice;
    
    cc = 0;
    for j = 1:size(ROI,1)
        ptest = ROI(j,:);
        [y, x] = meshgrid((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1));
        window = zeros(9,2);
        window(:,1) = y(:);
        window(:,2) = x(:);
        window(:,3) = slice; 
        
        window_intensity = post_sub((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1), slice);
        window_intensity = window_intensity(:);
        window_timecourse = mean(mean(DCE_LD_permuted((ptest(1)-1):(ptest(1)+1),(ptest(2)-1):(ptest(2)+1),:, slice)));
        window_timecourse = window_timecourse(:);

        cc_update = corr(window_intensity, kernel_intensity);
        cc_update = corr(window_timecourse, kernel_timecourse);        
        if cc_update > cc
            P = ptest;
            kernel = window;
            cc = cc_update;
        end
    end
    
    kernel_intensity = post_sub((P(1)-1):(P(1)+1), (P(2)-1):(P(2)+1), slice);
    kernel_intensity = kernel_intensity(:);
    kernel_timecourse = mean(mean(DCE_LD_permuted((P(1)-1):(P(1)+1),(P(2)-1):(P(2)+1),:, slice)));
    kernel_timecourse = kernel_timecourse(:);
    
    artery_right = [artery_right; kernel];
end
    
%%%%% down  
slice = slice0;
P = P0_right;
kernel_intensity = post_sub((P0_right(1)-1):(P0_right(1)+1), (P0_right(2)-1):(P0_right(2)+1), slice);
kernel_intensity = kernel_intensity(:);
kernel_timecourse = mean(mean(DCE_LD_permuted((P0_right(1)-1):(P0_right(1)+1),(P0_right(2)-1):(P0_right(2)+1),:, slice)));
kernel_timecourse = kernel_timecourse(:);

for i = 1:(slice0-1)
    slice = slice-1;
    [Y, X] = meshgrid((P(1)-2):(P(1)+2), (P(2)-2):(P(2)+2));
    ROI = zeros(25,3);
    ROI(:,1) = Y(:);
    ROI(:,2) = X(:);
    ROI(:,3) = slice;
    
    cc = 0;
    for j = 1:size(ROI,1)
        ptest = ROI(j,:);
        [y, x] = meshgrid((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1));
        window = zeros(9,2);
        window(:,1) = y(:);
        window(:,2) = x(:);
        window(:,3) = slice; 
        
        window_intensity = post_sub((ptest(1)-1):(ptest(1)+1), (ptest(2)-1):(ptest(2)+1), slice);
        window_intensity = window_intensity(:);
        window_timecourse = mean(mean(DCE_LD_permuted((ptest(1)-1):(ptest(1)+1),(ptest(2)-1):(ptest(2)+1),:, slice)));
        window_timecourse = window_timecourse(:);
        
        cc_update = corr(window_intensity, kernel_intensity);
        cc_update = corr(window_timecourse, kernel_timecourse);
        if cc_update > cc
            P = ptest;
            kernel = window;
            cc = cc_update;
        end
    end
    
    kernel_intensity = post_sub((P(1)-1):(P(1)+1), (P(2)-1):(P(2)+1), slice);
    kernel_intensity = kernel_intensity(:);
    kernel_timecourse = mean(mean(DCE_LD_permuted((P(1)-1):(P(1)+1),(P(2)-1):(P(2)+1),:, slice)));
    kernel_timecourse = kernel_timecourse(:);
%     
    artery_right = [artery_right; kernel];
end

%% remove unphysical voxels
remove = [];
artery = [artery_left;artery_right];
%
for i = 1:size(artery,1)
    timecourse = squeeze(DCE_LD_permuted(artery(i,1), artery(i,2), :, artery(i,3)));
    tc_3consecutive = conv2(timecourse, [1;1;1], 'valid')/3;
    [m, ind] = max(tc_3consecutive);
    
    pre_con = mean(timecourse(1:injtime)); %1->injection time 
    pre_std = std(timecourse(1:injtime)); %std(1->injection time)
    last_con = mean(timecourse((end-3):end)); %27->30
    
    enh_time = find((timecourse(6:end)-pre_con) > 5*pre_std); %end of injection time->30 > 5*pre_std
    if isempty(enh_time)  %remove all points that are not > 5*pre_std (enhancement)
         remove = [remove, i];
    else
        enh_time = enh_time(1)+5;
        post_con = mean(timecourse(enh_time:end));

        if  (m-pre_con)<5*pre_con|| (post_con-pre_con)<20*pre_std|| (last_con-pre_con)>0.85*(m-pre_con) ...   || (lastmax_con-pre_con) >0.85 * (m-pre_con)
             || pre_con <=0
            remove = [remove, i];
        else
            if  enh_time >enh0 || ind > enh_time+2
               remove = [remove, i];
            end
        end
        
    end
    
end

artery(remove,:) =[];
%}
%%% individual AIF SI
%%
if ~isempty(artery)
    timecourses = [];

    for i = 1:size(artery,1)
        timecourses = [timecourses, squeeze(DCE_LD_permuted(artery(i,1), artery(i,2), :, artery(i,3)))];
    end

    if ~isempty(timecourses)
        AIF = sum(timecourses,2) / size(artery,1);

        timecourses_SS0 = timecourses ./ repmat(mean(timecourses(1:5,:)), [numdyn,1]);
        AIF_intensity = sum(timecourses_SS0,2) / size(artery,1);
    end

%%%% plot results
%
    % figure(1)
    % plot(3.4*(1:length(AIF_intensity)), AIF_intensity)
    % xlabel('time(s)')
    % ylabel('SI/S_{pre}')
    % 
    % figure(2)
    % plot(3.4*(1:length(AIF_intensity)), timecourses_SS0)
    % xlabel('time(s)')
    % ylabel('SI/S_{pre}')
    figure(2)
    SIwithremoved = timecourses-repmat(mean(timecourses(1:5,:)),[numdyn,1]);
    numberofrows = size(SIwithremoved,1);
    numberofcolumns = size(SIwithremoved,2);
   
    AIFSI = [];
    AIFSI_index = [];
    for i = 1:numberofcolumns
        columntested = SIwithremoved(:,i);
        if isempty(find(columntested >= 1000000)) == false
            %SIwithremoved(:,i) = [];
            %SIwithremoved(:,i) = zeros(numberofrows,1);
            AIFSI = [AIFSI, SIwithremoved(:,i)];
            AIFSI_index = [AIFSI_index i];
        end
    end
  
    arterynew = zeros(91,3);
    for i = 1:length(AIFSI_index)
        arterynew(i,:) = (artery(AIFSI_index(i),:));
    end
    
    plot((1:length(AIF_intensity)), AIFSI)
    xlabel('time(s)')
    ylabel('SI-S_pre')

    AIFSI = AIFSI';
    editedAIF = mean(AIFSI);
    
    figure(3)
    plot((1:length(AIF_intensity)), editedAIF)
    xlabel('time(s)')
    ylabel('SI')

    figure(4)
    plot((1:length(AIF_intensity)), timecourses)
    xlabel('time(s)')
    ylabel('SI')    

end
%}    

%% Plotting different AIFS
AIFSItransposed = AIFSI'; %(time)x(# of SI time courses)
maxAIFSItransposed = max(AIFSItransposed); %(1)x(# of SI time courses)
% Determine number of values in:
% Top 10%:
top10percent = floor(0.1*size(AIFSItransposed,2));
% Top 15%:
top15percent = floor(0.15*size(AIFSItransposed,2));
% Top 20%:
top20percent = floor(0.2*size(AIFSItransposed,2));
% Top 25%:
top25percent = floor(0.25*size(AIFSItransposed,2));

% Sort peaks of SI time courses
[B,I] = sort(maxAIFSItransposed, 'descend'); %B = sorted array in descending order, I = original indices of maxAIFSItransposed

figure(5)
chosen_artery10 = [];
top10 = B([1:top10percent]); %array with top 10 maximum SI values
for i = 1:length(top10)
    chosen_artery10 = [chosen_artery10; arterynew(I(i),:)]; %choose SI time course
    chosentimecourse = AIFSItransposed(:,I(i)); %choose SI time course
    plot([1:1:size(AIFSItransposed,1)],chosentimecourse);
    hold on
    xlabel('time points');
    ylabel('SI(arb. units)');
    title('Top 10% AIFs');
end

figure(6)
chosen_artery15 = [];
top15 = B([1:top15percent]); %array with top 10 maximum SI values
for i = 1:length(top15)
    chosen_artery15 = [chosen_artery15; arterynew(I(i),:)]; %choose SI time course
    chosentimecourse = AIFSItransposed(:,I(i)); %choose SI time course
    plot([1:1:size(AIFSItransposed,1)],chosentimecourse);
    hold on
    xlabel('time points');
    ylabel('SI(arb. units)');
    title('Top 15% AIFs');
end

figure(7)
chosen_artery20 = [];
top20 = B([1:top20percent]); %array with top 10 maximum SI values
for i = 1:length(top20)
    chosen_artery20 = [chosen_artery20; arterynew(I(i),:)]; %choose SI time course
    chosentimecourse = AIFSItransposed(:,I(i)); %choose SI time course
    plot([1:1:size(AIFSItransposed,1)],chosentimecourse);
    hold on
    xlabel('time points');
    ylabel('SI(arb. units)');
    title('Top 20% AIFs');
end

figure(8)
chosen_artery25 = [];
top25 = B([1:top25percent]); %array with top 10 maximum SI values
for i = 1:length(top25)
    chosen_artery25 = [chosen_artery25; arterynew(I(i),:)]; %choose SI time course
    chosentimecourse = AIFSItransposed(:,I(i)); %choose SI time course
    plot([1:1:size(AIFSItransposed,1)],chosentimecourse);
    hold on
    xlabel('time points');
    ylabel('SI(arb. units)');
    title('Top 25% AIFs');
end

%% save data
%save([matpath subname '_Artery_tracking','.mat'], 'artery');

%% show image to check tracked artery
%{
test_image = post_sub;
for i = 1:size(artery,1)
    test_image(artery(i,1), artery(i,2), artery(i,3)) = 1000;
end

figure(2)
% colormap gray
for i = 10:(zsize-10)
    subplot(1,2,1)
    imagesc(test_image(:,:,i))
    caxis([0 800])
    axis equal tight
    title(num2str(i))
    
    subplot(1,2,2)
    imagesc(post_sub(:,:,i))
    caxis([0 800])
    axis equal tight
    pause(0.3)
    
    saveas(gcf, ['./', num2str(i), '.jpg'])
end
%}

%% save('Artery_right.mat', 'artery')
% 
% %% for both sides
% %
% load([matpath [subname '_Artery_left.mat']])
% artery_left = artery;
% 
% load([matpath [subname '_Artery_right.mat']])
% artery_right = artery;
% 
% artery = [artery_left;artery_right];

% %
% savepath = '/Users/chengyuewu/Desktop/study&work/Lab/VesselTracking/Workspace/Data/Processed/';
% save([savepath, subname, '_Artery_tracking.mat'],'artery', 'numdyn');
% %}


