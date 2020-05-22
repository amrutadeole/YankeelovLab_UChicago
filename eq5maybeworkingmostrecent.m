
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
%______________________
% Plot each pair of ve and ktrans values
for i = 1:length(assigned_ve)
    x = [assigned_ktrans(i), assigned_ve(i)];
    c_toi = kety_tofts(x,xd);
    %figure(i)
    %plot(t_new,c_toi, '.');
    c_toi_values = [c_toi_values c_toi];
end

%_________________________


%ktrans = (0.12/60);
%ve = 0.33;
%Ti = 0.40;

Cp = AIF_pop;
t = t_new;
R1i = (2/3); %s^-1
r1 = 4.5; %mM^-1*s^-1
R10 = (2/3); %s^-1
fw = 0.8;

x = [ktrans, ve];
xd = [t Cp];



for j = 1:10
    for i = 1:length(t_new)

    %R_t = 0.5*(2*R1i+(R10-R1i+1/Ti)/(ve/fw))-0.5*((2/Ti-(R10-R1i+1/Ti)/(ve/fw))^2+4*(1-ve/fw)/(Ti^2*(ve/fw)^0.5));
    %disp(R_t)
        A = (R10-R1i+(1/Ti(j)));
        B = ve(j)/fw;
        partA = 2*R1i+r1*c_toi(i,j)+A/B;
        partB = ((2/Ti(j))-r1*c_toi(i,j)-A/ve(j)/fw)^2;
        partC = 4*(1-B)/((Ti(j)^2)*(B^0.5));
        R(i) = (0.5*partA-0.5*(partB+partC));
        %R(i) = 0.5*(2*R1i+r1*c_toi(i)+A/B)-0.5*((2/Ti-r1*c_toi(i)-A/B)^2+4*((1-B)/Ti^2*B^0.5));
       % R(i) = 0.5*(2*R1i+r1*ktrans*crpexp_integral+(R10-R1i+1/assigned_Ti(i))/(assigned_ve(i)/fw))-0.5*((2/assigned_Ti(i)-r1*ktrans*crpexp_integral-(R10-R1i+1/assigned_Ti(i))/(assigned_ve(i)/fw))^2+4*(1-assigned_ve(i)/fw)/(assigned_Ti(i)^2*(assigned_ve(i)/fw)^0.5));
    end
end
R = R';
plot(t_new,R);

