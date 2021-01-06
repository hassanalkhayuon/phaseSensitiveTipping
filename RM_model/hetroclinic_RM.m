% the equlibria and stability of RM model.

% for detailes about the model and parameters go to:
% https://www.overleaf.com/project/5db3114c13676a00018d4902

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')

% close all
warning off
clear
clc

% parameters:

Rphas= 2.60865712;        %/yr.             Prey intrinsic growth rate
C     = 0.19;        %/yr.             Nonlinear death rate of prey
gamma = 0.004;       %pred/prey.       Prey-predator conversion rate
delta = 2.2;         %/yr.             Predator death rate in absent of prey
beta  = 1.5;         % prey/ha         Predator half-saturating constant
alpha = 800;         % prey/(pred.yr)  Predator saturating kill rate
mu    = 0.03;        % NA              Allee parameter
nu    = 0.003 ;      % NA              Allee parameter

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
%% phase portraits


odefun  = @(t,var)RModefun(var,Rphas);


e0_phas    = [0,0];
e1_phas    = [Rphas/C, 0];
e2_phas    = [mu,0];
e3_phas(1) = (delta*beta)/((gamma*alpha) - delta);
e3_phas(2) = (Rphas-(C*e3_phas(1)))*(beta+e3_phas(1))*(e3_phas(1)-mu)...
    /(alpha*(nu + e3_phas(1)));

% the stable manifold of e1
% maninit = e1_phas1 + [-0.01,0];
 maninit1 = e1_phas + [0 1e-8];
 maninit2 = e2_phas + [0 1e-8];
 maninit3 = e2_phas + [0.001 0];

[~,varman_phas1]   = ode45(odefun,[0,22.26],maninit1,opts);
[~,varman_phas2]   = ode45(odefun,[7.3,0],maninit2,opts);
[~,varman_phas3]   = ode45(odefun,[0,10],maninit3,opts);

figure(1);
cla
hold on
box on

plot(...
    e0_phas(1),e0_phas(2),'o',...
    'MarkerSize',3,'LineWidth',3.5,'Color','k')
plot(...
    e1_phas(1),e1_phas(2),'+',...
    e2_phas(1),e2_phas(2),'+',...
    'MarkerSize',7,'LineWidth',1.5,'Color','k')  
plot(...
    e3_phas(1),e3_phas(2),'o',...
    'MarkerSize',5,'LineWidth',2,'Color','k')

% plot(varper_phas(1,:),varper_phas(2,:),'-k','LineWidth',2)
% 
% plot(varper_phas(1,:),varper_phas(2,:),'Color',[.07 .6 1],'LineWidth',2)

plot(...
    varman_phas1(:,1),varman_phas1(:,2),'-b','LineWidth',2)
plot(...
    varman_phas2(:,1),varman_phas2(:,2),'-r','LineWidth',2)
plot(...
    varman_phas3(:,1),varman_phas3(:,2),'-r','LineWidth',2)


axis([-0.5 16 -.001 0.03])
set(gca,'FontSize',14)
xlabel('$R$')
ylabel('$P$','Rotation',0)
% 
%% functions

function [dvar] = RModefun(var,RR)
%parameters

C     = 0.19;        %/yr.             Nonlinear death rate of prey
gamma = 0.004;       %pred/prey.       Prey-predator conversion rate
delta = 2.2;         %/yr.             Predator death rate in absent of prey
beta  = 1.5;         % prey/ha         Predator half-saturating constant
alpha = 800;         % prey/(pred.yr)  Predator saturating kill rate
mu    = 0.03;        % NA              Allee parameter
nu    = 0.003 ;      % NA              Allee parameter

N = var(1);
P = var(2);

dN = RR * N *(1-((C/RR)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = gamma * ((alpha*P*N)/(beta + N)) - (delta*P);

dvar = [dN;dP];
end
