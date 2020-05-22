
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
assigned_ktrans = [0.3];
assigned_ktrans = assigned_ktrans./60;
assigned_ve = [0.15];
assigned_Ti = [1.0+2*(0.5/9)];

TR = 0.01;
flip = 5;

dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);


%R1_vox = R_toi_values;

% LSQ fit
x0st = [(0.1/60) 0.05 1.0];
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',1500,'MaxFunEval',2000);
for ii = 3
        ydata = S(:,ii);
        [x,resnorm] = lsqcurvefit('SIequation',x0st,t_new,ydata,[],[],options,Cp,R1i,r1,R10,fw,flip,TR,S0);
        calculated_ktrans(ii) = x(1)
        calculated_ve(ii) = x(2)
        calculated_Ti(ii) = x(3)
        %SIequation(x,t,Cp,flip,TR,R10,r1,fw,S0,R_toi_values)
        tempfit = SIequation([x(1) x(2) x(3)],t_new,Cp,R1i,r1,R10,fw,flip,TR,S0);
        fits(ii,:) = real(tempfit);
end

% Plot data and fit on graphs

    plot(t_new,S(:,3), '.');
    xlabel('Time(sec)');
    ylabel('SI(arb. units');
    hold on
    plot(t_new,fits(3,:));
    hold off



function [SI] = SIequation(x,t,Cp,R1i,r1,R10,fw,flip,TR,S0)
% Tofts equation and solve with
% linear least squared

%convert to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

n_images = numel(Cp);
Cp(isnan(Cp))= Cp(2);
R_1(1) = R10;
SI(1) = S0.*dcesin.*(1-exp(-TR.*R10))/(1-dcecos.*exp(-TR.*R10)); % SPGE eqn.

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
    
    %R_1(k)=0.5*(2*R1i+r1*c_toi(k)+(R10-R1i+(1/x(3)))/(x(2)/fw))-0.5*(((2/x(3))-r1*c_toi(k)-(R10-R1i+(1/x(3)))/(x(2)/fw))^2+4*(1-(x(2)/fw))/((x(3)^2)*(x(2)/fw)))^0.5;

    R_1(k) = 0.5*(2*R1i+r1*c_toi(k)+(R10-R1i+(1/(x(3))))/((x(2))/fw))-0.5*(((2/(x(3)))-r1*c_toi(k)-(R10-R1i+(1/(x(3))))/((x(2))/fw))^2+4*(1-((x(2))/fw))/(((x(3))^2)*((x(2))/fw)))^0.5;
    SI(k)= S0.*dcesin.*(1-exp(-TR.*R_1(k)))/(1-dcecos.*exp(-TR.*R_1(k)));


end

%SI = double(SI)
SI = SI';
%R1 = double(R1);
end