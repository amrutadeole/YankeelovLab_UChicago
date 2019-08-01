function [c_toi] = kety_tofts(x,xd)
% Tofts equation and solve with
% linear least squared

t = xd(:,1);
Cp = xd(:,2);

n_images = numel(Cp);
c_toi=zeros(n_images,1);
for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((x(1)/x(2).*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);%compute integral
    end
    t2 = t(1:k);%bounds for integral
    crpexp_integral = trapz(t2,crpexp);
    c_toi(k) = x(1)*crpexp_integral; %whole expression for R
  
end
c_toi = c_toi.';
end
