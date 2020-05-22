function S = t1_MFA_cf_findAIF(params,xd,alpha,TR, fitnoise,r1,T10mapvoxel)      
%Inputs:
%   -params: initial guesses [So ktrans ve]
%   -xd: AIF [t Cp]
%   -T10mapvoxel: T10 value for chosen voxel
%Outputs:
%   -signal intensity for each voxel
%   -fitted So, ktrans, ve parameters

%% Extract values from inputs
t = xd(:,1);
Cp = xd(:,2);
So = params(1);
ktrans = params(2);
ve = params(3);

%% Solve for the Ct(t) using Kety-Tofts Model
n_images = numel(Cp);
for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((ktrans/ve.*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);
    end
    t2 = t(1:k);%bounds for integral
    crpexp_integral = trapz(t2,crpexp); %compute integral
    c_toi(k,1) = ktrans*crpexp_integral; %whole expression for R
end

% Set v = 0 (there is no fitnoise)
if fitnoise == 1
    v = params(4);
else
    v = 0;
end

%% Solve for T1(t) term for each voxel--will fit for So, ktrans, ve
for i = 1:length(%SOME T TERM FOR AIF))
    T1_t_inv(t) = r1*c_toi(t)+1/T10mapvoxel; %T1mapvoxel is the chosen voxel's T10 value
    T1_t(t) = 1/T1_t_inv(t); 
    S(t) = double(So .* (sind(alpha).* (1 - exp(-TR/T1_t(t)))) ./ (1- (cosd(alpha).*exp(-TR/T1_t(t)))) +v);
 
end

end

    
    

