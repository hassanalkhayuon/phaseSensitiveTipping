% climate_tip_RM
% simulate one time series that tips to exctinction. 
% del is going to be fixes through out, and R is going to be varied btween
% Rstar and Rend.

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar      = 2.5;
Rend       = 1.6;
Rmean      = (Rstar+Rend)/2;
del        = 2.2;    %/yr.            Predator death rate in absent of prey
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
% rng(86) 
rng(1150) %for Rstar = 2.5
%%
Tend  = 100;
RR    = 1;
PnBin = 0.2;
% [swich,cli] = Climate(Tend,RR,PnBin);
[time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
climate = @(tt)interp1(time,cli1,tt);
R       = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;

figure
subplot(6,12,[1:7,13:19,25:31])
hold on
plot(time,Rmean + cli1.*((Rstar-Rend)/2),':k','LineWidth',1)
for ind = 1:length(swich)
    X = linspace(sum(swich(1:ind-1)),sum(swich(1:ind)),swich(ind));
    Y = Rmean + cli(sum(swich(1:ind-1))+1:sum(swich(1:ind)))*(Rstar-Rend)/2;
    if Y(1)>=Rmean
        plot(X,Y,'r-','LineWidth',2)
    else
        plot(X,Y,'b-','LineWidth',2)
    end
    if swich(ind) ==1
        X1 = [X-1,X];
        Y1 = [Y, Y];
        if Y1(1)>=Rmean
            plot(X1,Y1,'r-','LineWidth',2)
        else
            plot(X1,Y1,'b-','LineWidth',2)
        end
    end
end
box on
axis([0 Tend 1 3.1])
hA=gca;
set(gca,'XMinorTick','on','YMinorTick','on')
hA.XAxis.MinorTickValues = 0:5:100;
xticks([0 20 40 60 80 100])
xticklabels([])
hA.YAxis.MinorTickValues = 1:.1:3;
yticks([1.5 2 2.5 3 3.5 4])
yticklabels({1.5 2 2.5 3 3.5})
plot([0 200],[Rend Rend ],':k','LineWidth',1)
plot([0 200],[Rstar Rstar],':k','LineWidth',1)
set(gca,'FontSize',20)
ylabel('$R(t)$','Rotation',0)

%% time sieries
mult     = 1000;
initcond = [11,0.004];
odefun   = @(t,var)RModefun_cli(t,var,R);
[t,var]  = ode45(odefun,[0:0.1:Tend],initcond,opts);
subplot(6,12,[37:43,49:55,61:67])
hold on
plot(...
    t,var(:,1),'-k','LineWidth',2)
plot(...
    t,mult*var(:,2),'-','LineWidth',2,'Color',[.65 .65 .65])
box on
axis([0 Tend -2 30])
hA=gca;
set(gca,'XMinorTick','on','YMinorTick','on')
hA.XAxis.MinorTickValues = 0:5:100;
xticks([0 20 40 60 80 100])
xticklabels([0 20 40 60 80])
hA.YAxis.MinorTickValues = -2:1:45;
yticks([0 10 20 30])
yticklabels([0 10 20])
set(gca,'FontSize',20)
xlabel('$t$')
ylabel('$\tilde{P},N$','Rotation',90)
% var(end,1)
% crit(I) = var(end,1);
% end
%%
load('data/1par_bif.mat')

R_bif   = 0:0.001:4.5;

R_H     = 1.526;
R_h     = 2.609;
R_T1    = mu*C;
R_T2    = del*beta/(gamma * alpha - del)*C;
H_ind   = find(abs(R_bif-R_H)<0.0001);
h_ind   = find(abs(R_bif-R_h)<0.0001);
T1_ind  = find(abs(R_bif-R_T1)<0.0005);
T2_ind  = find(abs(R_bif-R_T2)<0.00025);


subplot(6,12,[8:12,20:24,32:36,44:48,56:60,68:72])
hold on

%1) R_1 R_2 lines 
plot([Rstar Rstar],[-10 20],':k','LineWidth',1)
plot([Rend Rend],  [-10 20],':k','LineWidth',1) 

% 2) equilibria curvies
e0      = zeros(2,length(R_bif));
e1(1,:) = R_bif./C;
e1(2,:) = zeros(1,length(R_bif));
e2(1,:) = mu*ones(1,length(R_bif));
e2(2,:) = zeros(1,length(R_bif));
e3(1,:) = (del*beta)/((gamma*alpha) - del) * ones(1,length(R_bif));
e3(2,:) = (R_bif -(C.*e3(1,:))).*(beta+e3(1,:)).*...
    (e3(1,:)-mu)./(alpha.*(nu + e3(1,:)));

plot(...
    R_bif,e0(1,:),'k-',...
    R_bif(T1_ind:T2_ind),e1(1,T1_ind:T2_ind),'k-',...
    R_bif(T2_ind:end),e1(1,T2_ind:end),'k--',...
    R_bif,e2(1,:),'k--',...
    R_bif(T2_ind:H_ind),e3(1,T2_ind:H_ind),'k-',...
    R_bif(H_ind:end),e3(1,H_ind:end),'k--',...
    'LineWidth',2)

%3) the periodic orbits 

D = 1;
n=1;
for ind_max =1:length(Rper)
    if rem(ind_max,D)==0
        Rpr(n)  = Rper(ind_max);
        Mxpr(n) = Maxper(ind_max);
        Mnpr(n) = Minper(ind_max);
        n=n+1;
    end
end

plot(Rpr,Mnpr,'-',Rpr,Mxpr,'-','LineWidth',2,'Color',[.47 .67 .19])

% %4) the periodic orbits 
plot(...
    R_bif(T1_ind),e1(1,T1_ind),'k.',...    % T1 point
    R_bif(T2_ind),e1(1,T2_ind),'k.',...  % T2 point
    R_bif(H_ind), e3(1,H_ind),'k.',...  % H point
    R_bif(h_ind), e1(1,h_ind),'k.',...  % H point
    'MarkerSize',30)

box on
axis([0 3 -1 20])
hA=gca;
set(gca,'XMinorTick','on','YMinorTick','on')
hA.XAxis.MinorTickValues = 0:.1:3;
xticks([0 1 2 3])
xticklabels([0 1 2 ])
hA.YAxis.MinorTickValues = 0:.5:20;
yticks([0 5 10 15 20])
yticklabels([0 5 10 15])
set(gca,'FontSize',20)
ylabel('$N$','Rotation',0)
xlabel('$R$')

% 5) zommed in 
axes('Position',[0.65 0.65 0.125 0.25])
hold on
box on 
set(gca,'FontSize',20)
xticklabels([])
yticklabels([])
plot(...
    R_bif,e0(1,:),'k-',...
    R_bif(1:T1_ind),e1(1,1:T1_ind),'k--',...
    R_bif(T1_ind:100),e1(1,T1_ind:100),'k-',...
    R_bif(1:T1_ind),e2(1,1:T1_ind),'-k',...
    R_bif(T1_ind:100),e2(1,T1_ind:100),'k--',...
    'LineWidth',2)
plot(...
    R_T1,e2(1,T1_ind),'k.',...    % T1 point
    'MarkerSize',30)
axis([0 0.01 0 0.053])

hA=gca;
set(gca,'XMinorTick','on','YMinorTick','on')
hA.XAxis.MinorTickValues = 0:.001:0.01;
xticks([0 0.005 0.01])
xticklabels([])
hA.YAxis.MinorTickValues = 0:0.005:0.055;
yticks([0 0.025 0.05])
yticklabels([])

%% functions
function [time,swich,climate,climate1] = Climate(Tend,RR,PnBin)
%climatswich generats a sequence of negative binomial random variables that
%add up to Tend, has propablility of success in singil tril PnBin and
% corresponding number of successes RR
swichtemp = nbininv(rand(1,2*Tend),RR,PnBin);
indend = 1;
while (sum(swichtemp(1:indend))<Tend)
    indend = indend+1;
end
swichtemp1 = swichtemp(1:indend);
ind_sw = 1;
for ind_sw1=1:length(swichtemp1)
    if swichtemp1(ind_sw1)~=0
        swich(ind_sw) = swichtemp1(ind_sw1);
        ind_sw = ind_sw+1;
    end
end
climate1 = NaN(1,sum(swich));
ind_cl = 1;
for ind_sw = 1:length(swich)
    climatampl = rand;
    if rand<=0.5 %mod(ind_sw,2)==1 % good years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = climatampl;
            ind_cl = ind_cl + 1;
        end
    else  % bad years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = -climatampl;
            ind_cl = ind_cl + 1;
        end
    end
end
time = 0:0.001:sum(swich);
climate = NaN(size(time));
for ind_cl=1:length(time)-1
    climate(ind_cl) = climate1(floor(time(ind_cl)+1));
end
end


%the ode fuction of RM model
function [dvar] = RModefun_cli(t,var,RRR)
%parameters

R          = RRR(t);
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


