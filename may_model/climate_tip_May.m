% climate_tip_May
% simulate one time series that tips to exctinction. 
% q is going to be fixes through out, and R is going to be decreased up
% to vary between Rstar and Rend and.

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar   = 3.3;       % /yr.              Prey intrinsic growth rate
Rend    = 2;
Rmean   = (Rstar+Rend)/2;
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-5,'AbsTol',1e-10);
%
rng(34)% rng(316) %rng(136)
%%
Tend  = 100;
RR    = 1;
PnBin = .2;
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
    if Y(1) >= Rmean
        plot(X,Y,'r-','LineWidth',2)
    else
        plot(X,Y,'b-','LineWidth',2)
    end
    if swich(ind) == 1
        X1 = [X-1,X];
        Y1 = [Y Y];
        if Y(1) >= Rmean
            plot(X1,Y1,'r-','LineWidth',2)
        else
            plot(X1,Y1,'b-','LineWidth',2)
        end
    end
end
box on
axis([0 Tend 1.5 4])
% hA=gca;
% set(gca,'XMinorTick','on','YMinorTick','on')
% hA.XAxis.MinorTickValues = 0:1:100;
% xticks([0 20 40 60 80 100])
% xticklabels([])
% hA.YAxis.MinorTickValues = 1.5:.1:4;
% yticks([1.5 2 2.5 3 3.5 4])
% yticklabels({'' 2 2.5 3 3.5})
% plot([0 150],[Rend Rend],':k','LineWidth',1)
% plot([0 150],[Rstar Rstar],':k','LineWidth',1)
% set(gca,'FontSize',15)
ylabel('$R(t)$','Rotation',0)
%% time sieries
mult     = 1000;
initcond = [11,0.004];
odefun  = @(t,var)Mayodefun_cli(t,var,R);
[t,var] = ode45(odefun,[0 Tend],initcond,opts);
subplot(6,12,[37:43,49:55,61:67])
hold on
plot(...
    t,var(:,1),'-k','LineWidth',2)
plot(...
    t,mult*var(:,2),'-','LineWidth',2,'Color',[.65 .65 .65])
box on
axis([0 Tend -2 38])
hA=gca;
% set(gca,'XMinorTick','on','YMinorTick','on')
% hA.XAxis.MinorTickValues = 0:5:100;
% xticks([0 20 40 60 80 100])
% xticklabels([0 20 40 60 80])
% hA.YAxis.MinorTickValues = -2:1:45;
% yticks([0 10 20 30])
% set(gca,'FontSize',15)
% xlabel('$t$')
% ylabel('$P,N$','Rotation',0)
%% bifurcation diagram
load('data/1par_bif.mat')
R_bif   = 0:0.005:4.5;

R_HB1 = 1.662;
R_HB2 = 3.81;
R_SN  = 1.205;
R_Fl  = 1.661;
HB_ind1 = find(abs(R_bif-R_HB1)<0.003);
HB_ind2 = find(abs(R_bif-R_HB2)<0.0001);
SN_ind  = find(abs(R_bif-R_SN)<0.001);

subplot(6,12,[8:12,20:24,32:36,44:48,56:60,68:72])
hold on

%1) R_1 R_2 lines 
% plot([Rstar Rstar],[-9 0],':k','LineWidth',1)
% plot([Rend Rend],[-9 0],':k','LineWidth',1) 

% 2) equilibria curvies 
plot(...
    R_bif,(e0(2,:)),'k-',...
    R_bif,(e1(2,:)),'k-',...
    R_bif,(e2(2,:)),'k-',...
    R_bif(1:HB_ind1),(e4(2,1:HB_ind1)),'k-',...
    R_bif(HB_ind1:end),(e4(2,HB_ind1:end)),'k--',...
    R_bif(1:HB_ind2),(e3(2,1:HB_ind2)),'k--',... 
    R_bif(HB_ind2:end),(e3(2,HB_ind2:end)),'k-',...
    R_bif,(e5(2,:)),'k--',...
    'LineWidth',2)


%3) the periodic orbits 

D = 1;

n=1;
for ind_max =1:length(Rper)
    if rem(ind_max,D)==0
        Rpr(n)  = Rper(ind_max);
        Mxpr(n) = Maxperp(ind_max);
        Mnpr(n) = Minperp(ind_max);
        n=n+1;
    end
end

% Min_per = @(RR)interp1(Rpr,Mnpr,RR);
% Max_per = @(RR)interp1(Rpr,Mxpr,RR);

% Rper = linspace(1.6618,3.8032,500);

plot(Rpr,(Mnpr),'-',Rpr,(Mxpr),'-','LineWidth',2,'Color',[.47 .67 .19])

%4) the periodic orbits 
plot(...
    R_bif(SN_ind),(e4(2,SN_ind)),'k.',...    % SN point
    R_bif(HB_ind1),(e4(2,HB_ind1)),'k.',...  % first HB point
    R_bif(HB_ind2),(e3(2,HB_ind2)),'k.',...  % second HB point
    'MarkerSize',25)

box on
% axis([1 4.05 -4 -.5])
hA=gca;
set(gca,'XMinorTick','on','YMinorTick','on')
hA.XAxis.MinorTickValues = 1:.1:4;
xticks([1 2 3 4])
xticklabels([1 2 3])
hA.YAxis.MinorTickValues = -9:.2:0;
yticks([-4 -3 -2 -1 0])
yticklabels([-4 -3 -2])
set(gca,'FontSize',15)
ylabel('$\log(\tilde{P})$','Rotation',90)
xlabel('$R$')

%%
% changing climate
function [time,swich,climate,climate1] = Climate(Tend,RR,PnBin)
%climatswich generats a sequence of negative binomial random variables that
%add up to Tend, has propablility of success in singil tril PnBin and
% corresponding number of successes RR
swichtemp = nbininv(rand(1,2*Tend),RR,PnBin); % how many swiches
indend = 1;
while (sum(swichtemp(1:indend))<Tend)
    indend = indend+1;
end
swichtemp1 = swichtemp(1:indend); % number of swiches befor Tend
ind_sw = 1;
for ind_sw1=1:length(swichtemp1)% clearing the swich list from zeros
    if swichtemp1(ind_sw1)~=0
        swich(ind_sw) = swichtemp1(ind_sw1);
        ind_sw = ind_sw+1;
    end
end
climate1 = NaN(1,sum(swich)); %climat vector sum(swich) long
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


% ode of the May Model
function [dvar] = Mayodefun_cli(t,var,RRR)

%parameters
R       = RRR(t);
q       = 205;
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.


N = var(1);
P = var(2);



dN = R * N *(1-((C/R)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = s*P*( 1-((q*P)/(N+eps)) );

dvar = [dN;dP];

end



