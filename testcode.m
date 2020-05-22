
% Clear workspace

% Load data
load('C:/Yankeelov Lab/AIF_Pop_TXO.mat')

% Initialize ve and ktrans values
assigned_ve = [0.23];
assigned_ktrans = [0.063];
t_new = t_new./60;
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


ktrans = (0.063);
ve = 0.23;
Ti = 0.91;

Cp = AIF_pop;
%t = t_new;
R1i = (2/3); %s^-1
r1 = 4.5; %mM^-1*s^-1
R10 = (2/3); %s^-1
fw = 0.8;

x = [ktrans, ve];
xd = [t_new Cp];


for i = 1:length(t_new)
    R_toi(i) = 0.5*(2*R1i+r1*c_toi_column_chosen(i)+(R10-R1i+(1/Ti))/(ve/fw))+0.5*(((2/Ti)-r1*c_toi_column_chosen(i)-(R10-R1i+(1/Ti))/(ve*fw))^2+4*(1-(ve/fw))/((Ti^2)*(ve/fw)))^0.5;
end
R = R';
plot(t_new,R_toi);

TR = 0.01;
flip = 5;

dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

for i = 1:10
    for j = 1:length(t_new)
        S(j,i) = S0.*dcesin.*(1-exp(-TR.*R1_vox(j,i)))/(1-dcecos.*exp(-TR.*R1_vox(j,i)));
    end
end

for i = 1:10
    plot(t_new,S(:,i))
end


