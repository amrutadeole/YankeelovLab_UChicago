function [R_toi] = equation5function(x,t,Cp,R1i,r1,R10,fw)
% Tofts equation and solve with
% linear least squared
% Create empty values to store each set of c(t) values 

load('savedvariables','c_toi_values')
%t_new = t;
%AIF_pop = Cp;
ktrans = x(1);
ve = x(2);
Ti = x(3);

n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);

R_toi = zeros(1,length(t));
R_toi(1) = R10;

for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((x(1)/x(2).*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);
    end
    t2 = t(1:k);
    crpexp_integral = trapz(t2,crpexp);
    R_toi(k) = 0.5*(2*R1i+r1*ktrans*crpexp_integral+(R10-R1i+1/Ti)/(ve/fw))-0.5*((2/Ti-r1*ktrans*crpexp_integral-(R10-R1i+1/Ti)/(ve/fw))^2+4*(1-ve/fw)/(Ti^2*(ve/fw)^0.5));
end


end