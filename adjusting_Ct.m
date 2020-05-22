
for i = 1:12
    Cp_for_adjusted(i) = 4.9618.*baseline_AIF(i);
end

for i = 13:86
    Cp_for_adjusted(i) = 318.5345*exp(-0.4720*i)+0.1727*exp(0.0076*i);
end

time_adjusted = [1:1:86];

%for i = 1:length(param_out)
i = 896145;
    concentration_adjusted(i,:) = kety_tofts(param_out(i,:)',[time_adjusted,Cp_for_adjusted]);
%end

