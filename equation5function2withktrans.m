function [R_toi] = equation5function2withktrans(x,t,Cp,R1i,r1,R10,fw)
% Tofts equation and solve with
% linear least squared
t_new = t;
AIF_pop = Cp;
ktrans = x(1);
ve = x(2);
Ti = x(3);


n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);

for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((ktrans/ve.*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);
    end
    t2 = t(1:k);
    crpexp_integral = trapz(t2,crpexp);
%     crpexp_integral = trapzfm(t2,crpexp); % quicker?
    c_toi_column_chosen(k) = x(1)*crpexp_integral;
end


%R_toi(1) = (2/3);
for i = 1:length(t_new)
    R_toi(i) = 0.5*(2*R1i+r1*c_toi_column_chosen(i)+(R10-R1i+(1/Ti))/(ve/fw))-0.5*(((2/Ti)-r1*c_toi_column_chosen(i)-(R10-R1i+(1/Ti))/(ve/fw))^2+4*(1-(ve/fw))/((Ti^2)*(ve/fw)))^0.5;
end


end