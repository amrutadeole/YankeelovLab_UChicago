function [SI] = SIequationwithnoise_FXL_7_29(x,t,Cp,flip,TR,R1map_sec,r1,s_inh,voxelindex)
% Tofts equation and solve with
% linear least squared

R10 = R1map_sec(voxelindex(1),voxelindex(2),voxelindex(3));
s_inh_voxel = s_inh(voxelindex(1),voxelindex(2),voxelindex(3));
% Convert to radians
sina = sind(flip);
cosa = cosd(flip);

n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);

% Initial values of R1 and SI
R1(1) = R10;
SI(1) = s_inh_voxel*sina*(1-exp(-TR*R10))/(1-cosa*exp(-TR*R10)); % SPGE eqn.

% Calculate Ct(t), R1(t), SI(t) values
for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((x(1)/x(2).*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);
    end
    t2 = t(1:k);
    crpexp_integral = trapz(t2,crpexp);
    c_toi(k) = x(1)*crpexp_integral;
    
    R1(k)=(r1*c_toi(k))+R10; 
    
    SI(k)= s_inh_voxel*sina*(1-exp(-TR*R1(k)))/(1-cosa*exp(-TR*R1(k)));
     %needs to be SI(k) = SI(k)+alpha*storerandom(k)

end
disp(c_toi)
SI = double(SI);
SI = SI';


end