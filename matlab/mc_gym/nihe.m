% Assume you determined xdata and ydata experimentally
% xdata = -6:0.1:10;
% xdata=xdata';
% 
% ydata = 100.*exp(-(xdata-2).^2./0.003^2)+10;
% noise = 4.8 .* rand(1, length(ydata));
% noise = noise';
% ydata = ydata + noise;
% % plot(xdata, ydata, 'o');
% % hold on;
% 
% fun = fittype('A*exp(-((x-u)/sigma)^2)+n');
% % fun = fittype('(A/(w*sqrt(PI/2)))*exp(-2*((x-xc)/w)^2)+(A1/(w1*sqrt(PI/2)))*exp(-2*((x-xc1)/w1)^2)');
% coeffnames(fun);  
% options = fitoptions();
% % options.StartPoint = [4 10 0.05 2];  % ?????
% 
% f=fit(xdata,ydata,'gauss2',options);
% fit1 = f(xdata);
% plot (xdata, fit1)

%% MATLAB
function lol = nihe(args)
    arg1 = args.arg1;
    arg2 = args.arg2;
    lol = arg1 + arg2;
end