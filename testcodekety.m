
% Tofts equation and solve with
% linear least squared

t = t_new
Cp = AIF_pop
x = [21/20000 0.23];

n_images = numel(Cp);

for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((x(1)/x(2).*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);%compute integral
    end
    t2 = t(1:k);%bounds for integral
    crpexp_integral = trapz(t2,crpexp);
    c_toi(k,1) = x(1)*crpexp_integral; %whole expression for R
  
end

