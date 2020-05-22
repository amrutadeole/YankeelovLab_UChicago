
   
%% Initialize values
TR = 4.820; %repetition time in ms
r1 = 4.5;
R1 = 1./(T1map); %in 1/seconds
R1(T1map(:)==0) = 0;
s_t = double(reshaped_DCE_LD_array);
flip = 30;

% Inherent Signal (S0 or M0), second 5 reps, before injection ~ index 12 (2
% min)
s_inh = mean(s_t(:,:,:,6:10),4).*(1-exp(-TR*R1).*cosd(flip))./((1-exp(-TR*R1))*sind(flip));
s_inh(isinf(s_inh)) = 0;
s_inh(isnan(s_inh)) = 0;

%% Convert 91x30 SI array into 91x30 R(t) array
for i = 1:91
    for t = 1:30
        s_inh_matching = s_inh(artery(AIFSI_index(i),1),artery(AIFSI_index(i),2),artery(AIFSI_index(i),3));
        R(i,t) = abs((-1/TR)*log((AIFSI(i,t)-s_inh_matching*sind(flip))./(AIFSI(i,t)*cosd(flip)-s_inh_matching*sind(flip))));
    end
end

%% Convert 91x30 R(t) array into 91x30 C(t) array
Ct = zeros(91,30);
for i = 1:91
    for t = 1:30
        T10matching = T1map(artery(AIFSI_index(i),1),artery(AIFSI_index(i),2),artery(AIFSI_index(i),3)); %choose corresponding T1 value  v
        Ct(i,t) = ((1/r1)*(R(i,t)-T10matching));
    end
end

for i = 24
    plot(1:length(AIF_intensity), Ct(i,:))
    hold on
    xlabel('time(sec)')
    ylabel('C(mM)')
end

%{
Ctrows = size(Ct,1);
Ctcolumns = size(Ct,2);
for i = 1:Ctrows
    for j = 1:Ctcolumns
        if Ct(i,j+1) - Ct(i,j) > 0.01
            fprintf("You have liftoff at %2.4f seconds\n", Ct(j+1));
        end
    end
end

%}