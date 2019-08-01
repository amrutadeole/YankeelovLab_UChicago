function [SI_t_FXL] = LDtermfunction(params,timerange_FXL,Cp_FXL,voxelindex,r1,flip,TR,s_inh_LD,R1map_sec)
    
    ktrans = params(1);
    ve = params(2);
    s_inh_voxel = s_inh_LD(voxelindex(1),voxelindex(2),voxelindex(3));
    R10 = R1map_sec(voxelindex(1),voxelindex(2),voxelindex(3));
    %R10 = R1map_sec(160,160,12);
    sina = sind(flip);
    cosa = cosd(flip);
    
    c_toi = kety_tofts([ktrans ve], [timerange_FXL' Cp_FXL']);
    for i = 1:length(c_toi)
        R1_t_FXL(i) = r1*c_toi(i)+R10; %in mM and sec
        
        SI_t_FXL(i) = s_inh_voxel*(sina*(1-exp(-TR*R1_t_FXL(i)))/(1-cosa*exp(-TR*R1_t_FXL(i))));
    end
end