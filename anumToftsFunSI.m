function [SI] = anumToftsFunSI(x,t,Cp,flip,TR,R1_vox,Rel,S0)
% Tofts equation and solve with
% linear least squared

%convert to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);
R1(1) = R1_vox;
SI(1) = S0.*dcesin.*(1-exp(-TR.*R1_vox))/(1-dcecos.*exp(-TR.*R1_vox)); % SPGE eqn.

for k = 2:n_images
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((x(1)/x(2).*(int_t-dummy_t))));
        crpexp(j) = Cp(j)*expo(j);
    end
    t2 = t(1:k);
    crpexp_integral = trapz(t2,crpexp);
%     crpexp_integral = trapzfm(t2,crpexp); % quicker?
    c_toi(k) = x(1)*crpexp_integral;
    R1(k)=(Rel*c_toi(k))+R1_vox; 
    %     R1(k)=Rel*(1/0.8)*c_toi(k)+R1_vox; % Some people set fraction of
    %     water accessible to CA to 0.8, some people set it to 1. Anum set it
    %     to 0.8 for preclinical Vandy data analysis, but will set it to 1 for
    %     all future UT analysis (clinical and preclinical)
    %if exist('noise','var')
     %   SI(k)= (1+random('Uniform',-1*noise,noise))*S0.*dcesin.*(1-exp(-TR.*R1(k)))/(1-dcecos.*exp(-TR.*R1(k)));
        SI(k)= S0.*dcesin.*(1-exp(-TR.*R1(k)))/(1-dcecos.*exp(-TR.*R1(k)));


end
SI = double(SI);
%R1 = double(R1);
end