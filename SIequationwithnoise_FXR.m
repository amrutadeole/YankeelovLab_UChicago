function [SI] = SIequationwithnoise_FXR(x,t,Cp,R1i,r1,R10,fw,flip,TR,S0,storerandom)
% Tofts equation and solve with
% linear least squared

%convert to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);
%R_1(1) = R10;
%SI(1) = S0.*dcesin.*(1-exp(-TR.*R10))/(1-dcecos.*exp(-TR.*R10)); % SPGE eqn.

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
    
end

for i = 1:length(t)
    %R_1(k)=0.5*(2*R1i+r1*c_toi(k)+(R10-R1i+(1/x(3)))/(x(2)/fw))-0.5*(((2/x(3))-r1*c_toi(k)-(R10-R1i+(1/x(3)))/(x(2)/fw))^2+4*(1-(x(2)/fw))/((x(3)^2)*(x(2)/fw)))^0.5;
    R_1(i) = 0.5*(2*R1i+r1*c_toi(i)+(R10-R1i+(1/(x(3))))/((x(2))/fw))-0.5*(((2/(x(3)))-r1*c_toi(i)-(R10-R1i+(1/(x(3))))/((x(2))/fw))^2+4*(1-((x(2))/fw))/(((x(3))^2)*((x(2))/fw)))^0.5;
    SI(i)= S0.*dcesin.*(1-exp(-TR.*R_1(i)))/(1-dcecos.*exp(-TR.*R_1(i)));
    SI(i) = SI(i)+0.26175*storerandom(i);
end

%SI = double(SI)
%SI = SI';
%R1 = double(R1);
end