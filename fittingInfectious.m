clc;clearvars;close all;
load('KCDC.mat');

close all
plot(T.Date,T.InfectedCum,'x')
hold on
plot(T.Date,T.Isolated,'x')
plot(T.Date,T.Discharged,'x')
legend('Cummlative Confirmed Total','Confirmed and Isolated','Confirmed and Dischaged') 

% f = fittype('rat33')

sdf(gcf,'paper_f150')

%%
% clc;clearvars
% % mybeta = 0.00038;
% % mygamma = 0.1;
% syms I(t) mybeta S I mygamma I;
% ode = diff(I,t) == mybeta*S*I-mygamma*I;
% iSol(t) = dsolve(ode)

%%
% % function didt = f(t,i)
% % didt = 
% 
% % https://kr.mathworks.com/matlabcentral/answers/476516-how-to-fit-differential-equations-to-a-curve
% syms c(t) k1 k2 t c0
% Eqn = diff(c) == 6*k1 - k1*t - k2*t^2;
% C = dsolve(Eqn, c(0) == c0)
% t = [5 10 20 30 45 60 90 120];
% c = [4.83 3.87 2.54 2.08 1.82 1.8 1.76 1.74];
% B = [ones(size(t(:))), 6*t(:)-t(:).^2/2, -t(:).^3/3] \ c(:)
% cf = [ones(size(t(:))), 6*t(:)-t(:).^2/2, -t(:).^3/3] * B;
% figure
% plot(t, c, 'p')
% hold on
% plot(t, cf, '-r')
% hold off
% grid
%% https://kr.mathworks.com/help/optim/ug/fit-differential-equation-ode.html
clc;close all
% sigma = 10;
% beta = 8/3;
% rho = 28;
% mybeta = 0.0000000018;
% mygamma = 0.01;
% population = 50000000;

mybeta = 0.00000078;
mygamma = 0.01;
population = 50000;

f = @(t,a) [-mybeta*a(1)*a(2); mybeta*a(1)*a(2) - mygamma*a(2); mygamma*a(2)];
xt0 = [population,1,0];
[tspan,a] = ode45(f,[1 230],xt0);     % Runge-Kutta 4th/5th order ODE solver
figure
% plot(a(:,2))
plot(a)
% legend('Susceptible ', 'Infected ', 'Recovered ');
% view([-10.0 -2.0])
%%
soln = T.Isolated;
soln(isnan(soln))=0;

xt0 = zeros(1,3);
xt0(1) = mybeta;
xt0(2) = mygamma;
xt0(3) = population;

[pbest,presnorm,presidual,exitflag,output] = lsqcurvefit(@paramfun,xt0,tspan,soln(floor(tspan)));
%%
fprintf('New parameters: %f, %f, %f\n',pbest(1:3))
fprintf('Original parameters: %f, %f, %f\n',xt0)

odesl = presidual + soln(floor(tspan));
plot(floor(tspan),odesl,'b')
hold on
plot(floor(tspan),soln(floor(tspan)))


function pos = paramfun(x,tspan)

mybeta = x(1);
mygamma = x(2);
f = @(t,a) [-mybeta*a(1)*a(2); mybeta*a(1)*a(2) - mygamma*a(2); mygamma*a(2)];
xt0 = [x(3),1,0];
[~,pos] = ode45(f,tspan,xt0);
pos = pos(:,2);
end
