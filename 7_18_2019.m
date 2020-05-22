
tic
%% Data Pre-processing
% Set no bounds
lb = [];
ub = [];
alpha = 30;
% Data
MFA_T1W_data = cat(4,FA3,FA5,FA10,FA15,FA20,FA30); %create 4D: (y,x,slice #,FA #)
MFA_T1W_data2 = MFA_T1W_data(:,:,12,:); %test just one slice
[m, n, o, p] = size(MFA_T1W_data2); %o is number of slices, p is number of FA scans
FAmap_ratio_data = ones([m,n,o]); %set FAmap_ratio_data = 1;
FAs = [3,5,10,15,20,30]; %array of FAs in ascending order

% Repetition time in ms    
TR = 4.820; 

%% Parameters for Fitting

im_data = reshape(MFA_T1W_data2, [m*n*o, p]); %reshaped data so that every voxel in 3D space is a row, and every FA is a column
FAmap_vec = reshape(FAmap_ratio_data, [m*n*o, 1]); %all equal to 1

% Calculate initial guesses
T1_init = mean(FAs);
So=max(MFA_T1W_data2(:));
params = [So 1000]; %params(1) = So*;  params(2) = tisue T1* guess;

param_out = zeros(m*n*o,length(params));

resid_map = zeros(m*n*o,1); %initialize map of residuals
fails = zeros(m*n*o,1); %initialize matrix for fails
options = optimset('TolFun',1e-6,'Tolx',1e-6,'MaxIter',500,'Display','off');%,'algorithm','levenberg-marquardt');

% Variables to improve efficiency
nodes = m*n*o;
parloops = nodes/(m*0.5); % # chunks at a time

%% Voxel wise calculation of T1 Map, parallelized
for jstep = 1:(nodes/parloops) % looping through chunks
    j_list = 1+parloops*(jstep-1):parloops*jstep; %list of nodes within a chunk
    
    im_data_pf = im_data(j_list,:);%chunk of im_data
    FAmap_vec_pf = FAmap_vec(j_list,:);
    param_out_pf = zeros(parloops,length(params));
    resid_map_pf = zeros(parloops,1);
    fails_pf = zeros(parloops,1);
    
    parfor ii = 1:parloops
        if ~all(im_data_pf(ii,:)==0) %&& FAmap_vec_pf(ii,:)~=0 % zero values for FA map or T1W data ignore
            im_vox = im_data_pf(ii, :);
            FAmap_vox = FAmap_vec_pf(ii,:); %ratio value
            FAs_scaled = FAmap_vox.*FAs;
         
            % Fitting
            try
                [param_out_pf(ii,:), resid_map_pf(ii,1)] = lsqcurvefit('t1_MFA_cf',params,FAs_scaled,im_vox,lb,ub,options,TR,0);
            catch
                fails_pf(ii,1) =1;
            end
        end
    end
    
    param_out(j_list,:) = param_out_pf;
    resid_map(j_list,:) = resid_map_pf;
    fails(j_list,:)= fails_pf;
    if sum(fails(:))>0
        warning('Warning: Had some fails!')
    end
    
    disp(['T1map ' num2str(round((jstep/(nodes/parloops))*100,2)) '% finished...'])
end


T1map = reshape(param_out(:,2),[m,n,o]);
Somap = reshape(param_out(:,1),[m,n,o]);

MFA_T1W_residual = reshape(resid_map(:,1),[m,n,o]);
fails_map = reshape(fails,[m,n,o]);

toc
