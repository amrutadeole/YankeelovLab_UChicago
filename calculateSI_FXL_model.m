% Forward Fitting Curve


assigned_ve = [0.05:0.05:0.5];
assigned_ktrans = [0.1:0.1:1.0];
assigned_ktrans = assigned_ktrans./60;

S0 = 364.45; %chosen value from S0 map
r1 = 4.5; %mM^-1*s^-1
R10 = (1/1.5); %s^-1
flip = 5; %degrees
TR = 0.01; %s

rows = size(c_toi_values,1); 
columns = size(c_toi_values,2);

% Convert flip angle to radians
dcesin = sin((flip/180)*pi);
dcecos = cos((flip/180)*pi);

% Initialize R1_vox matrix
R1_vox = zeros(rows,columns);

% R1_vox matrix for each ktrans and ve combination
for i = 1:columns
   for j = 1:rows
       R1_vox(j,i) = r1*c_toi_values(j,i)+R10;
   end
end

% Solve for signal intensity for each R value
for i = 1:columns
    for j = 1:rows
        S(j,i) = S0.*dcesin.*(1-exp(-TR.*R1_vox(j,i)))/(1-dcecos.*exp(-TR.*R1_vox(j,i)));
    end
end


% LSQ fit

x0st=[0.1 0.2];
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',150,'MaxFunEval',2000);
for ii = 1:size(S,2)        
        ydata = S(:,ii);
        [x,resnorm] = lsqcurvefit('amrutaToftsFunSI',x0st,t_new,ydata,[],[],options,AIF_pop,flip,TR,R10,r1,S0);
        calculated_ktrans(ii) = x(1);
        calculated_ve(ii) = x(2);
        tempfit = amrutaToftsFunSI([x(1) x(2)],t_new,AIF_pop,flip,TR,R10,r1,S0);
        fits(ii,:) = tempfit;
end


% Plot data and fit on graphs
for jj = 1:10
    figure(jj)
    plot(t_new,S(:,jj), '.');
    xlabel('Time(sec)');
    ylabel('SI(arb. units');
    str = sprintf('SI Vs. Time Curve with ktrans = %4.4f and ve = %4.2f', assigned_ktrans(jj), assigned_ve(jj)); title(str);
    hold on
    plot(t_new,fits(jj,:));
    ylim([20,35]);
    hold off
end

% Plot family of curves on same graph
for jj = 1:10
    figure(11)
    plot(t_new,S(:,jj));
    xlabel('Time(sec)');
    ylabel('SI(arb. units)');
    title('Family of Curves');
    hold on
end

 
 


