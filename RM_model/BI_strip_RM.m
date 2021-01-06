
warning off
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');



% parameters:

Rlow       = 2;   %/yr.            Prey intrinsic growth rate
Rupp       = 2.5;
Rend       = Rupp;
del        = 2.2;    %/yr.            Predator death rate in absent of prey
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter

T        = 100;
 
 opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
%% step 1 defining periodic orbits

mult     = 1000;
initcond = [8   0.01];
Tper     = 10;%8.1480;

ivpfun   = @(t,var)RModefun(var,Rlow);
perfun   = @(t,var,Tper)Tper*RModefun(var,Rlow);

[~,tempvar]  = ode45(ivpfun,[0 T],initcond);
initper  = tempvar(end,:); %the top point is [6.9261, 0.0339];
ivpsol = ode45(@(t,var)perfun(t,var,Tper),[0 1],initper, opts);


% The boundary condtions for the periodic orbit
ss = linspace(0,1,200);
tempinit = @(s)deval(ivpsol,s);
solinit=bvpinit(ss,tempinit,Tper);

BC=@(WL,WR,Tper)(...
    [WL(1)-WR(1);...
    WL(2)-WR(2);...
    WL(2)-initper(2);... %point phase condition
    ]);

per_low = bvp5c(perfun,BC,solinit);
load('data\perinitstart_RM1.mat')

per_upp = BVPper_start;

N_initconds = 1000; % number of initial conditions on the limit cycle
ss = linspace(0,1,N_initconds);

varlow = deval(per_low,ss);
varupp = deval(per_upp,ss);


figure; hold on
plot(...
    varlow(1,:)+1.2,mult*varlow(2,:)+3.3,...
    'LineWidth',3,'Color',[0.47 0.67 0.19])
plot(...
    varupp(1,:),mult*varupp(2,:),...
    'LineWidth',2,'Color',[0.47 0.67 0.19])
%%
Rlow       = 2; 
ivpfun   = @(t,var)RModefun(var,Rlow);
[~,var]  = ode45(ivpfun,[0 9],[10.5 0.0001],opts);
% [~,var]  = ode45(ivpfun,[0 12],[3.31    8.1557/1000+0.001],opts);
figure
plot(...
    var(:,1),mult*var(:,2),...
    'LineWidth',1.5,'Color','k')
%% replicating the two periodic orbits and \theta(R)

% the start periodic orbit
Nupp_top = [(varupp(1,506:1000)),(varupp(1,2:104))];
Pupp_top = mult*[(varupp(2,506:1000)),(varupp(2,2:104))];

Nupp_dow = varupp(1,104:506);
Pupp_dow = mult*varupp(2,104:506);

Pstar_top = @(NN)interp1(Nupp_top,Pupp_top,NN);
Pstar_dow = @(NN)interp1(Nupp_dow,Pupp_dow,NN);


% the end periodic orbit

Nlow_top = varlow(1,119:632);
Plow_top = mult*varlow(2,119:632);

Nlow_dow = [varlow(1,632:999),varlow(1,1:119)];
Plow_dow = mult*[varlow(2,632:999),varlow(2,1:119)];

Pend_top = @(NN)interp1(Nlow_top,Plow_top,NN);
Pend_dow = @(NN)interp1(Nlow_dow,Plow_dow,NN);

% 2 threshold 
odefun  = @(t,var)RModefun(var,Rend);
e2      = [mu,0];
G       = @(var)RModefun(var,Rend);
JJ      = MyJacobian(G,e2);
[eigvic,eigval] = eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-4*pert';
maninit = e2 + pert;
[~,varman]   = ode45(odefun,[10,0],maninit,opts);
Pthe = @(NN)interp1(varman(:,1),mult*varman(:,2),NN);

%% scan in (N,P)-plane

grid_N = 1000;
grid_P = 1.5*grid_N;
N_scan = linspace(0,15,grid_N);
P_scan = linspace(0,25,grid_P);
colormat = zeros(grid_P,grid_N);

for ind_N = 1:grid_N
    N  = N_scan(ind_N);
        
%     % Region 2.1
        Pmax2_1 = min(Pthe(N),Pstar_top(N));
        Pmin2_1 = Pend_top(N);
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_1, P>Pmin2_1)
                colormat(ind_P,ind_N) = 2;
            end
        end
     
    % Region 2.2
    if or(N>=5.195,N<=1.922)
        Pmax2_2 = Pstar_top(N);
        Pmin2_2 = Pstar_dow(N);
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_2, P>Pmin2_2)
                colormat(ind_P,ind_N) = 2;
            end
        end 
    end
    
    % Region 2.3
    if or(N<=5.195,N>=1.922)
        Pmax2_3 = Pend_dow(N);
        Pmin2_3 = Pstar_dow(N);
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_3, P>Pmin2_3)
                colormat(ind_P,ind_N) = 2;
            end
        end 
    end
    % Region 1

        Pmax1 = Pstar_top(N);
        Pmin1 = max(Pthe(N),Pend_top(N));
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax1, P>Pmin1)
                colormat(ind_P,ind_N) = 1;
            end
        end
end

%% plotting 
load('data\freq_phase_RM_v2_rand.mat')
figure;
map = [1,1,1;1 0 0;.47,.67,.19];
PCOLOR = pcolor(N_scan,P_scan,colormat);
PCOLOR.FaceAlpha = 0.1;
PCOLOR.LineStyle = 'none';
colormap(map)
hold on
plot(var_bef(1,:),mult*var_bef(2,:),'k.','Markersize',7)
% fplot(Pstar_top,[0 15])
% fplot(Pstar_dow,[0 15])
% fplot(Pend_top,[0 15])
% fplot(Pend_dow,[0 15])
fplot(Pthe,[0,15],'--k','LineWidth',1)
caxis([0 2]);

%%  arrows 
% Xi =[7.4 6.62 6.351];
% Yi = [17.55 17.6 17.13]
% X= [min(Xi):0.01:max(Xi)];
% Y = interp1(Xi,Yi,X,'spline');
% plot(X,Y,'LineWidth',1.1,'Color',[.65 .65 .65])

%% functions 

%the ode fuction of RM model
function [dvar] = RModefun(var,par)
%parameters
% global C gamma beta alpha mu nu delta

RR    = par(1);


delta      = 2.2;
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter

N = var(1);
P = var(2);

dN = RR * N *(1-((C/RR)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = gamma * ((alpha*P*N)/(beta + N)) - (delta*P);

dvar = [dN;dP];
end

