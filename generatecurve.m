% Clear workspace

% Load data
load('C:/Yankeelov Lab/AIF_Pop_TXO.mat')

% Initialize ve and ktrans values
assigned_ve = [0.05:0.05:0.5];
assigned_ktrans = [0.1:0.1:1.0];
assigned_ktrans = assigned_ktrans./60;
t_new = t_new';
AIF_pop = AIF_pop';
xd = [t_new AIF_pop];

% Create empty values to store each set of c(t) values 
c_toi_values = [];

% Plot each pair of ve and ktrans values
for i = 1:length(assigned_ve)
    x = [assigned_ktrans(i), assigned_ve(i)];
    c_toi = kety_tofts(x,xd);
    %figure(i)
    %plot(t_new,c_toi, '.');
    c_toi_values = [c_toi_values c_toi];
end

assigned_ktrans = assigned_ktrans';
assigned_ve = assigned_ve';

% Initial guess
x0=[0.1, 0.05];
    
% LSQ fit
for ii = 1:size(c_toi_values,2)        
        ydata = c_toi_values(:,ii);
        [x,resnorm] = lsqcurvefit(@kety_tofts,x0,xd,ydata);
        calculated_ktrans(ii) = x(1);
        calculated_ve(ii) = x(2);
        tempfit = kety_tofts([x(1) x(2)], xd);
        fits(ii,:) = tempfit;
end

% Plot data and fit on graphs
for jj = 1:10
    figure(jj)
    plot(t_new,c_toi_values(:,jj), '.');
    xlabel('Time(sec)');
    ylabel('Ct(mM)');
    str = sprintf('Ct Vs. Time Curve with ktrans = %4.4f and ve = %4.2f', assigned_ktrans(jj), assigned_ve(jj)); title(str);
    hold on
    plot(t_new,fits(jj,:));
    ylim([0,2]);
    hold off
end
    
% Plot family of curves on same graph
for jj = 1:10
    figure(11)
    plot(t_new,c_toi_values(:,jj));
    xlabel('Time(sec)');
    ylabel('Ct(mM)');
    title('Family of Curves');
    hold on
end

calculated_ktrans = calculated_ktrans';
calculated_ve = calculated_ve';
assigned_ktrans = assigned_ktrans.*60;
calculated_ktrans = calculated_ktrans.*60;

% Create table with ktrans and ve values
T = table(assigned_ktrans, assigned_ve, calculated_ktrans, calculated_ve)
