ve = 0.1;
ktrans = 0.139/60; %delete later
%TR = 4.820/1000;
t = [1:1:30];
x0st=[0.2];
options = optimset('TolFun',1e-12,'TolX',1e-12);options = optimset(options,'LargeScale','on');
options = optimset(options,'Display','off');options = optimset(options,'MaxIter',150,'MaxFunEval',2000);
        
ydata = AIF_gluteal_2';
%[x, resnorm] = lsqcurvefit('scalingAIF_prostatefunction',x0st,t,ydata,[],[],options,AIF_prostate_tenpercent, ve);
[x, resnorm] = lsqcurvefit('scalingAIF_prostatefunction',x0st,t,ydata,[],[],options,AIF_prostate_tenpercent,ktrans,ve);
%calculated_ktrans = x(1);
%calculated_alpha = x(2);
calculated_alpha = x(1)
