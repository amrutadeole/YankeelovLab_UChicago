       

% Load data
load('C:/Yankeelov Lab/AIF_Pop_TXO.mat')
load('eq5variables')


Cp = AIF_pop;
t = t_new;
R1i = (2/3); %s^-1
r1 = 4.5; %mM^-1*s^-1
R10 = (2/3); %s^-1
fw = 0.8; %fraction value

% Initial values
assigned_ktrans = [0.1:0.1:1.0];
assigned_ktrans = assigned_ktrans./60;
assigned_ve = [0.05:0.05:0.5];
assigned_Ti = [1.0:(0.5/9):1.5];


% x is a 3x1 column vector with [Ktrans; ve; Ti]
n_images = numel(Cp); %count the number of elements in AIF_pop
%Cp(isnan(Cp))= Cp(2);

R_toi_values = [];

for i = 1:length(assigned_ve)
    x = [assigned_ktrans(i) assigned_ve(i) assigned_Ti(i)];
    c_toi_column_chosen = c_toi_values(:,i);
    R_toi = equation5function2(x,t,Cp,R1i,r1,R10,fw,c_toi_column_chosen);
    R_toi = R_toi';
    R_toi_values = [R_toi_values R_toi];
    plot(t_new,R_toi_values(:,i));
end

format long
% LSQ fit
x0st = [(0.1/60) 0.05 1.0];
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',1500,'MaxFunEval',2000);
%temp_fit_values = [];
for ii = 1:size(R_toi_values,2) %for each column in R_toi_values       
        ydata = R_toi_values(:,ii);
        ydata = ydata';
        c_toi_column_chosen = c_toi_values(:,ii);
        [x,resnorm] = lsqcurvefit(@equation5function2,x0st,t,ydata,[],[],options,Cp,R1i,r1,R10,fw,c_toi_column_chosen);
        calculated_ktrans(ii) = x(1);
        calculated_ve(ii) = x(2);
        calculated_Ti(ii) = x(3);
        tempfit = equation5function2([x(1) x(2) x(3)],t_new,AIF_pop,R1i,r1,R10,fw,c_toi_column_chosen);
        %tempfit = tempfit';
        %temp_fit_values = [temp_fit_values tempfit];
        fits(ii,:) = tempfit;
end

% Plot data and fit on graphs
for jj = 1:10
    figure(jj)
    plot(t,R_toi_values(:,jj), '.');
    ylim([0.6,1.4]);
    xlabel('Time(sec)');
    ylabel('R(t) (1/sec)');
    str = sprintf('R(t) Vs. Time Curve with ktrans = %4.4f, ve = %4.2f, and Ti = %4.2f', assigned_ktrans(jj), assigned_ve(jj), assigned_Ti(jj)); title(str);
    hold on
    plot(t,fits(jj,:));
    hold off
end

% Plot family of curves on same graph
for jj = 1:10
    figure(11)
    plot(t_new,R_toi_values(:,jj));
    xlabel('Time(sec)');
    ylabel('Rt(t) (1/sec)');
    title('Family of Curves');
    hold on
end

assigned_ktrans = real(assigned_ktrans)';
assigned_ve = real(assigned_ve)';
assigned_Ti = real(assigned_Ti)';
calculated_ktrans = real(calculated_ktrans)';
calculated_ve = real(calculated_ve)';
calculated_Ti = real(calculated_Ti)';
T = table(assigned_ktrans, assigned_ve, assigned_Ti, calculated_ktrans, calculated_ve, calculated_Ti)