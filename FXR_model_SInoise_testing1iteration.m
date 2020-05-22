% Outputs:
%   - SI_3d = 3-dimensional array with:
%       - (t=66) rows
%       - (runs = 100 for each combo) columns
%       - (combos of ktrans, ve, and Ti = 10) slices 
%   - calculated_ktrans =  100x10 matrix 
%       - (runs = 100 for each combo) rows
%       - (combos of ktrans, ve, and Ti = 10) columns

% Load data

tic

load('C:/Yankeelov Lab/AIF_Pop_TXO.mat')
load('SIvariables')

% Define variables
t = t_new;
R1i = (2/3); %s^-1
r1 = 4.5; %mM^-1*s^-1
R10 = (2/3); %s^-1
fw = 0.8; %fraction value
S0 = 364.45; %chosen value from S0 map
TR = 0.01;
flip = 5;

% Initial values
assigned_ktrans = 0.1;
assigned_ktrans = assigned_ktrans/60;
assigned_ve = 0.05;
assigned_Ti = 1.0;

storerandom = normrnd(0,1);


% Calcualate radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

% Initialize arrays 
SI_values = [];
SI_3d = zeros(66,1,1);  

% Forward fitting 
 % i = "vertical slice" number
    x = [assigned_ktrans assigned_ve assigned_Ti]; %a new vertical "slice" of the cube
    % k = # of runs per "vertical slice"
        SI = SIequationwithnoise_FXR(x,t,Cp,R1i,r1,R10,fw,flip,TR,S0, storerandom);  
        SI = SI';
        SI_3d(:,1,1) = SI;

% LSQ fit
x0st = [(0.4/60) 0.35 1.0];
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',4000,'MaxFunEval',4000);
fits = zeros(66,1,1);


        ydata = SI_3d(:,1,1);
        ydata = ydata';
        [x,resnorm] = lsqcurvefit('SIequationwithktrans',x0st,t_new,ydata,[],[],options,Cp,R1i,r1,R10,fw,flip,TR,S0);
        calculated_ktrans = x(1);
        calculated_ve = x(2);
        calculated_Ti = x(3);
        tempfit = SIequationwithktrans([x(1) x(2) x(3)],t_new,Cp,R1i,r1,R10,fw,flip,TR,S0);
        fits(:,1,1) = real(tempfit);


%{
for i = 1:10
    mean_ktrans(i) = sum(calculated_ktrans(:,i))./100;
    mean_ve(i) = sum(calculated_ve(:,i))./100;
    mean_Ti(i) = sum(calculated_Ti(:,i))./100;
end

assigned_ktrans = assigned_ktrans';
assigned_ve = assigned_ve';
assigned_Ti = assigned_Ti';
mean_ktrans = mean_ktrans';
mean_ve = mean_ve';
mean_Ti = mean_Ti';    

format short
%T = table(assigned_ktrans, assigned_ve, calculated_ktrans, calculated_ve)
Table = table(assigned_ktrans, assigned_ve, assigned_Ti, mean_ktrans, mean_ve, mean_Ti)

%calculate standard deviation of each parameter for each combination
for i = 1:10
    vector = calculated_ktrans(:,i);
    std_ktrans(i) = std(vector);
    vector = calculated_ve(:,i);
    std_ve(i) = std(vector);
    vector = calculated_Ti(:,i);
    std_Ti(i) = std(vector);
end

for i = 1:10
    alpha_ktrans(i) = std_ktrans(i)/10;
    alpha_ve(i) = std_ve(i)/10;
    alpha_Ti(i) = std_Ti(i)/10;
end
    
%calculate means of each parameter
mean_ktrans = sum(calculated_ktrans)./100;
mean_ve = sum(calculated_ve)./100;
mean_Ti = sum(calculated_Ti)./100;

figure(1) 
plot(t_new, SI_3d(:,1,1),'.');
hold on
plot(t_new, fits(:,1,1));
toc

% NEED TO DO PARALLEL PROCESSING FOR GRAPHING
for iii = 1:10
    figure (iii)
    parfor jjj = 1:100
        loop_ktrans = [0.1:0.1:1.0];
        loop_ktrans = assigned_ktrans./60;
        loop_ve = [0.05:0.05:0.5];
        loop_Ti = [1.0:(0.5/9):1.5];
        subplot(10,10,jjj)
        plot(t_new, SI_3d(:,jjj,iii), '.'); %plot generated data points
        ylabel('SI(arb. units)');
        xlabel('Time(sec)');
        str = sprintf('SI Vs. Time Curve with ktrans = %4.4f, ve = %4.2f, and Ti = %4.2f', loop_ktrans(iii), loop_ve(iii), loop_Ti(iii)); title(str);
    end
end

%toc


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
%}
