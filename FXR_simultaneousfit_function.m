function [SI] = SIequationwithnoise_FXR_7_30(params,timerange_FXR,Cp_FXR,r1,fw,flip,TR,s_inh,R1map_sec,voxelindex)
% Tofts equation and solve with
% linear least squared
%load('storerandom')
%convert to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

s_inh_voxel = s_inh(voxelindex(1),voxelindex(2),voxelindex(3));

t = timerange_FXR;
n_images = numel(Cp_FXR);
Cp_FXR(isnan(Cp_FXR))= Cp_FXR(2);
%R_1(1) = R10;
%SI(1) = S0.*dcesin.*(1-exp(-TR.*R10))/(1-dcecos.*exp(-TR.*R10)); % SPGE eqn.
R10 = R1map_sec(voxelindex(1),voxelindex(2),voxelindex(3));
R1i = R10;

for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((params(1)/params(2).*(int_t-dummy_t))));
        crpexp(j) = Cp_FXR(j)*expo(j);
    end
    t2 = t(1:k);
    crpexp_integral = trapz(t2,crpexp);
    c_toi(k) = params(1)*crpexp_integral;
    
end

for i = 1:length(t)
    R_1(i) = 0.5*(2*R1i+r1*c_toi(i)+(R10-R1i+(1/(params(3))))/((params(2))/fw))-0.5*(((2/(params(3)))-r1*c_toi(i)-(R10-R1i+(1/(params(3))))/((params(2))/fw))^2+4*(1-((params(2))/fw))/(((params(3))^2)*((params(2))/fw)))^0.5;
    SI(i)= s_inh_voxel.*dcesin.*(1-exp(-TR.*R_1(i)))/(1-dcecos.*exp(-TR.*R_1(i)));
    %SI(i) = SI(i)+0.26175*storerandom(i);
end

%SI = double(SI)
%SI = SI';
%R1 = double(R1);
end
