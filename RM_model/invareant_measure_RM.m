% frequency of tipping in terms of phases \varphi
warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rinv   = 2.2;       %/yr.            Prey intrinsic growth rate
del    = 2.2;       %/yr.            Predator death rate in absent of prey
C      = 0.19;      %/yr.            Nonlinear death rate of prey
gamma  = 0.004;     %pred/prey.      Prey-predator conversion rate
beta   = 1.5;       % prey/ha        Predator half-saturating constant
alpha  = 800;       % prey/(pred.yr) Predator saturating kill rate
mu     = 0.03;      % NA             Allee parameter
nu     = 0.003;     % NA             Allee parameter

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);

%% scan (N,P)-plane

einv(1) = (del*beta)/((gamma*alpha) - del);
einv(2) = (Rinv-(C*einv(1)))*(beta+einv(1))*(einv(1)-mu)/(alpha*(nu + einv(1)));

mult = 1000;
odefun   = @(t,var)RModefun(var,Rinv);

K = 10000;

Nscan = linspace(0,20,K);
Pscan = linspace(0,30,K);
% figure;
% hold on
% endPoint = NaN(2,K^2);
indTotal = 1;
for indN = 1:K
    for indP = 1: K
        initCondition = [Nscan(indN),Pscan(indP)/mult];
        [t,var] = ode45(odefun,[0,50],initCondition,opts);
        if norm(var(end,:))>1e-3
            endPoint(:,indTotal) = var(end,:);
            indTotal = indTotal + 1;
        end
    end
    disp([indN,K])
end

[phiinv,rinv]  = cart2pol(endPoint(1,:)-einv(1),...
    (mult.*(endPoint(2,:)-einv(2))));
figure
histogram(...
    phiinv,60,'Normalization','count',...
    'FaceColor','[0.47 0.67 0.19]'...
    ,'LineStyle','-','LineWidth',1.5)

%%

%the ode fuction of RM model
function [dvar] = RModefun(var,R)
%parameters

delta      = 2.2;
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter


N = var(1);
P = var(2);

dN = R * N *(1-((C/R)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = gamma * ((alpha*P*N)/(beta + N)) - (delta*P);

dvar = [dN;dP];
end
