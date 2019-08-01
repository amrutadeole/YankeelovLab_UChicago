function [ct_muscle] = scalingAIF_prostatefunction(x,t,Cp,ktrans,ve)

% Code for scaling AIF(artery) to AIF(muscle) in prostate
%   -fits ktrans and alpha values
%   -x = [ktrans_initialguess, alpha_initialguess]
%   -Cp = AIF from artery
%   -t = time point range

AIF_prostate_new = Cp;
n_images = numel(AIF_prostate_new);
for k = 2:n_images %2:30
    int_t = t(k);
    for j = 1:k
        dummy_t = t(j);
        expo(j) =exp(-((ktrans/ve.*(int_t-dummy_t))));
        crpexp(j) = x(1)*AIF_prostate_new(j)*expo(j);
    end
    t2 = t(1:k);%bounds for integral
    crpexp_integral = trapz(t2,crpexp); %compute integral
    ct_muscle(k,1) = ktrans*crpexp_integral; %whole expression for R
end

