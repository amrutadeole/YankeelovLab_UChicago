
% Load data
load('C:/Yankeelov Lab/AIF_Pop_TXO.mat')
load('SIvariables')

t = t_new;
R1i = (2/3); %s^-1
r1 = 4.5; %mM^-1*s^-1
R10 = (2/3); %s^-1
fw = 0.8; %fraction value
S0 = 364.45; %chosen value from S0 map

% Initial values
assigned_ktrans = [0.1:0.1:1.0];
assigned_ktrans = assigned_ktrans./60;
assigned_ve = [0.05:0.05:0.5];
assigned_Ti = [1.0:(0.5/9):1.5];

TR = 0.01;
flip = 5;

dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);


SI_values = [];

for i = 1:length(assigned_ve)
    x = [assigned_ktrans(i) assigned_ve(i) assigned_Ti(i)];
    %c_toi_column_chosen = c_toi_values(:,i);
    SI = SIequationwithktrans(x,t,Cp,R1i,r1,R10,fw,flip,TR,S0);
    SI = SI';
    SI_values = [SI_values SI];
end


% LSQ fit
x0st = [(0.4/60) 0.35 1.2];
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',4000,'MaxFunEval',4000);
for ii = 1:10
        ydata = SI_values(:,ii);
        ydata = ydata';
        [x,resnorm] = lsqcurvefit('SIequationwithktrans',x0st,t_new,ydata,[],[],options,Cp,R1i,r1,R10,fw,flip,TR,S0);
        calculated_ktrans(ii) = x(1);
        calculated_ve(ii) = x(2);
        calculated_Ti(ii) = x(3);
        tempfit = SIequationwithktrans([x(1) x(2) x(3)],t_new,Cp,R1i,r1,R10,fw,flip,TR,S0);
        fits(ii,:) = real(tempfit);
end



% Plot data and fit on graphs
for jj = 1:10
    figure(jj)
    plot(t_new,SI_values(:,jj), '.');
    ylim([20,25])
    xlabel('Time(sec)');
    ylabel('SI(arb. units)');
    str = sprintf('SI Vs. Time Curve with ktrans = %4.4f, ve = %4.2f, and Ti = %4.2f', assigned_ktrans(jj), assigned_ve(jj), assigned_Ti(jj)); title(str);
    hold on
    plot(t_new,fits(jj,:));
    hold off
end

% Plot family of curves on same graph
for jj = 1:10
    figure(11)
    plot(t_new,SI_values(:,jj));
    xlabel('Time(sec)');
    ylabel('SI(arb. units)');
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