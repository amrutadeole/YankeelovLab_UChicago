function [R_vox] = calculatervox(t,Cp,x,index,assigned_ktrans,assigned_ve,assigned_Ti)

    R1i = (1/1.5); %s^-1
    r1 = 4.5; %mM^-1*s^-1
    R10 = (1/1.5); %s^-1
    fw = 0.8;
    clear crpexp
    R_vox(1) = (1/1.5);
    
    for p = 2:t
        %k = length(t);
        %n_images = numel(Cp);

        for k = 2:p
            int_t = t(k);
            for j = 1:k
                dummy_t = t(j);
                expo(j) =exp(-((x(1)/x(2).*(int_t-dummy_t))));
                crpexp(j) = Cp(j)*expo(j);%compute integral
            end
            
            t2 = t(1:k);%bounds for integral
            crpexp_integral = trapz(t2,crpexp);
            R_vox(k) = 0.5*(2*R1i+r1*assigned_ktrans(index)*crpexp_integral+(R10-R1i+1/assigned_Ti(index)/(assigned_ve(index)/fw))-0.5*(2/assigned_Ti(index)-r1*assigned_ktrans(index)*crpexp_integral-(R10-R1i+1/assigned_Ti(index))/(assigned_ve(index)/fw)^2)+4*(1-assigned_ve(index)/fw)/(assigned_Ti(index)^2*(assigned_ve(index)/fw)^0.5));
        end
    end
end