
   
%% Initialize values
TR = 4.820; %repetition time in ms
r1 = 4.5/1000;
R1 = 1./(T1map); %in 1/ms
R1(T1map(:)==0) = 0;
s_t = double(reshaped_DCE_LD_array);
flip = 30;

% Inherent Signal (S0 or M0), second 5 reps, before injection ~ index 12 (2
% min)
s_inh = mean(s_t(:,:,:,6:10),4).*(1-exp(-TR*R1).*cosd(flip))./((1-exp(-TR*R1))*sind(flip));
s_inh(isinf(s_inh)) = 0;
s_inh(isnan(s_inh)) = 0;

%% Convert 91x30 SI array into 91x30 R(t) array
for i = 1:432
    for t = 1:30
        s_inh_matching = s_inh(artery(i,1),artery(i,2),artery(i,3));
        R_lessfilter(i,t) = abs((-1/TR)*log((plotAIF(i,t)-s_inh_matching*sind(flip))./(plotAIF(i,t)*cosd(flip)-s_inh_matching*sind(flip))));
    end
end

%% Convert 91x30 R(t) array into 91x30 C(t) array
Ct_lessfilter = zeros(432,30);
for i = 1:432
    for t = 1:30
        T10matching = T1map(artery(i,1),artery(i,2),artery(i,3)); %choose corresponding T1 value  v
        Ct_lessfilter(i,t) = ((1/r1)*(R_lessfilter(i,t)-T10matching));
    end
end

for i = 1:432
    plot(1:length(AIF_intensity), Ct_lessfilter(i,:))
    hold on
    xlabel('time(sec)')
    ylabel('C(mM)')
end

%}
