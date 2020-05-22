
load('C:/Yankeelov Lab/AIF_Pop_TXO.mat')
load('allvariables')

r1 = 4.5; %mM^-1*s^-1
R10 = (1/1.5); %s^-1
R1i = (1/1.5); %s^-1
t = t_new;
Cp = AIF_pop;
fw = 0.8;

assigned_ktrans = [0.1:0.1:1.0];
assigned_ve = [0.05:0.05:0.5];
assigned_Ti = [1.5:(0.25/9):1.75];

n_images = numel(Cp);

Rt_t_values = [];
for i = 1:size(c_toi_values,2)
    c_toi_vector = c_toi_values(:,i);
    Rt_t = calculateRt_t(x,t,Cp,R1i,r1,c_toi_vector,fw,R10);
    Rt_t_values = [Rt_t_values Rt_t];
    figure(i);
    plot(t_new,Rt_t);
end

function [Rt_t] = calculateRt_t(x,t,Cp,R1i,r1,c_toi_vector,fw,R10)

    temp_ktrans = x(1);
    temp_ve = x(2);
    temp_Ti = x(3);
    Rt_t(1) = (2/3);
    
    for i = 2:length(t)
        Rt_t(i) = 0.5*(2*R1i+r1*c_toi_vector(i)+(R10-R1i+(1/temp_Ti))/(temp_ve/fw))-0.5*(((2/temp_Ti)-r1*c_toi_vector(i)-(R10-R1i+(1/temp_Ti))/temp_ve/fw)^2+4*(1-(temp_ve/fw))/(temp_Ti^2)*(temp_ve/fw)^0.5);
    end
    

end


