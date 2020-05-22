function [SI] = amrutaToftsFunSI(x,t,Cp,flip,TR,R10,Rel,S0)
% Tofts equation and solve with linear least squared


%convert to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);

% Initial R1 and SI values
R1(1) = R10;
SI(1) = S0.*dcesin.*(1-exp(-TR.*R10))/(1-dcecos.*exp(-TR.*R10)); % SPGE eqn.

% Calculate Ct(t), R1(t), and SI(t) values
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
    
    R1(k)=(Rel*c_toi(k))+R10; 

    SI(k)= S0.*dcesin.*(1-exp(-TR.*R1(k)))/(1-dcecos.*exp(-TR.*R1(k)));

end
SI = double(SI);
SI = SI';
%R1 = double(R1);
end