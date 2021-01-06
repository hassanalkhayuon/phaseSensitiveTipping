% the equlibria and stability of RM model.

% for detailes about the model and parameters go to:
% https://www.overleaf.com/project/5db3114c13676a00018d4902

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')

% close all
warning off
clear
clc

% parameters:
global C gamma beta alpha mu nu

R     = 0:0.001:3.5;  %/yr.             Prey intrinsic growth rate
C     = 0.19;        %/yr.             Nonlinear death rate of prey
gamma = 0.004;       %pred/prey.       Prey-predator conversion rate
delta = 2.2;         %/yr.             Predator death rate in absent of prey
beta  = 1.5;         % prey/ha         Predator half-saturating constant
alpha = 800;         % prey/(pred.yr)  Predator saturating kill rate
mu    = 0.03;        % NA              Allee parameter
nu    = 0.003 ;      % NA              Allee parameter

% In Rebecca simulations mu = 0.45 and  and nu = 0.15;

opts = odeset('RelTol',1e-4,'AbsTol',1e-8);

e0 = NaN(2,length(R));
e1 = NaN(2,length(R));
e2 = NaN(2,length(R));
e3 = NaN(2,length(R));

for nn = 1:length(R)
    RR = R(nn);
    G=@(var)RModefun(var,RR);
    
    e0(1,nn) = 0;
    e0(2,nn) = 0;
    
    e1(1,nn) = RR/C;
    e1(2,nn) = 0;
    
    e2(1,nn) = mu;
    e2(2,nn) = 0;
    
    e3(1,nn) = (delta*beta)/((gamma*alpha) - delta);
    N         = e3(1,nn);
    e3(2,nn) = (RR-(C*N))*(beta+N)*(N-mu)/(alpha*(nu + N));
    P         = e3(2,nn);
    
    JJ = MyJacobian(G,e0(:,nn));
    if isfinite(sum(sum(JJ)))
        stability0(:,nn)=eig(JJ);
    else
        stability0(:,nn) = NaN(2,1);
    end
    
    JJ = MyJacobian(G,e1(:,nn));
    if isfinite(sum(sum(JJ)))
        stability1(:,nn)=eig(JJ);
    else
        stability1(:,nn) = NaN(2,1);
    end
    
    JJ = MyJacobian(G,e2(:,nn));
    if isfinite(sum(sum(JJ)))
        stability2(:,nn)=eig(JJ);
    else
        stability2(:,nn) = NaN(2,1);
    end
    
    JJ = MyJacobian(G,e3(:,nn));
    if isfinite(sum(sum(JJ)))
        stability3(:,nn)=eig(JJ);
    else
        stability3(:,nn) = NaN(2,1);
    end
    
end

%  transcritical points

R_tra1 = mu*C;
R_tra2 = delta*beta/(gamma * alpha - delta)*C;

N_tra1 = mu;
N_tra2 = R_tra2/C;

P_tra1=0;
P_tra2=0;

% hopf point

R_H = 1.53;

N_H = (delta*beta)/((gamma*alpha) - delta);
P_HB = (R_H-(C*N_H))*(beta+N_H)*(N_H-mu)/(alpha*(nu + N_H));

ind_tra1 = find(abs(R - R_tra1)<0.01);
ind_tra2 = find(abs(R - R_tra2)<0.01);
ind_HB   = find(abs(R-R_H)<0.01);

% homo point

R_h = 2.609;
N_h = R_h/C;
P_h = 0;

%% 2 parameter diagram

load('RM2par_new.mat')

del  = linspace(0.3765,3.1,500);
del1 = linspace(0,3.1,500);
R_tra1_2par  = C*mu*(ones(size(del1)));
R_tra2_2par = C*(del1*beta)./(gamma*alpha - del1);
% load('2pardia_xppaut.mat') % the homoclinic and the Hopf from XPPaut.
TT=[C*mu,mu*alpha*gamma/(mu+beta)];
GH=[0.2741 0.5673];
Rhind = 219;
Rh = h2par(Rhind,:);
%% plotting
mult = 1000;
figure
set(0,'defaulttextInterpreter','latex')
hold on

%1) the bifurcation points

plot3(R_tra1,N_tra1,mult*P_tra1,'k.',...
    R_tra2,N_tra2,mult*P_tra2,'k.',...
    R_H,N_H,mult*P_HB,'k.',...
    R_h,N_h,mult*P_h,'k.',...
    'MarkerSize',20);

% 2) orgin equilibrium (always stable)
plot3(R,e0(1,:),mult*e0(2,:),'k','LineWidth', 1.5)

% 3)e1 stable between R=R_tra1 and R= R_tra2;
plot3(R(ind_tra1(1):ind_tra2(1)),e1(1,ind_tra1(1):ind_tra2(1)),...
    mult*e1(2,ind_tra1(1):ind_tra2(1)),'-k','LineWidth',1.5)
plot3(R(ind_tra2(1):end),e1(1,ind_tra2(1):end),...
    mult*e1(2,ind_tra2(1):end),'--k','LineWidth',1.5)

% 4) e2, stable between  R=0 and R=R_tra1
plot3(R(1:ind_tra1(1)),e2(1,1:ind_tra1(1)),mult*e2(2,1:ind_tra1(1)),...
    '-k','LineWidth',1.5)
plot3(R(ind_tra1(1):end),e2(1,ind_tra1(1):end),mult*e2(2,ind_tra1(1):end),...
    '--k','LineWidth',1.5)

% 5) e3, stabel between R = R_tra2 and R=R_HB
plot3(R(ind_tra2:ind_HB),e3(1,ind_tra2:ind_HB),mult*e3(2,ind_tra2:ind_HB),...
    '-k','LineWidth',1.5)
plot3(R(ind_HB:end),e3(1,ind_HB:end),mult*e3(2,ind_HB:end),...
    'k--','LineWidth',1.5)

% 6) plotting the periodic orbits
load('RM_per.mat')
ss=linspace(0,1,500);
for kk = 1:10
    varper = deval(per(kk),ss);
    plot3(Rper(kk)*ones(1,500),varper(1,:),mult*varper(2,:),'-k')
end

axis([0 3 0 18 0 mult*0.035])
xlabel('$R$')
ylabel('$N$')
zlabel('$P$','Rotation',0)
view([-70,40])
set(gca,'FontSize',14)
grid on
box on

annotation('textbox',...
    [0.3367 0.496 0.0782 0.092],...
    'String','$h$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14);

annotation('textbox',...
    [0.658 0.315 0.076 0.0773],...
    'String','$H$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','on');
annotation('textbox',...
    [0.5117 0.183 0.0957 0.0923],...
    'String','$T_2$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','on');
annotation('textbox',...
    [0.2688 0.2489 0.0784 0.0771],...
    'String','$T_1$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','on');
set(gca,'FontSize',14)

% 
% annotation('arrow',[0.6679 0.3357],...
%     [0.1276 0.2786],...
%     'Color',[0.35 0.75 0.9]);

% zoomed in

R_zoom = [linspace(0,R_tra1,5),linspace(R_tra1,0.015,5)];

for mm = 1:length(R_zoom)
    RR = R_zoom(mm);
    
    e0_zoom(1,mm) = 0;
    e0_zoom(2,mm) = 0;
    
    e1_zoom(1,mm) = RR/C;
    e1_zoom(2,mm) = 0;
    
    e2_zoom(1,mm) = mu;
    e2_zoom(2,mm) = 0;
    
end

%%
axes('Position',[0.13 0.255 0.25 0.25])
box on
hold on
%1) the bifurcation points
plot3(R_tra1,N_tra1,mult*P_tra1,'k.',...
    'MarkerSize',24);

% 2) e0 equilibrium (always stable)

plot3(R_zoom,e0_zoom(1,:),mult*e0_zoom(2,:),'-k','LineWidth', 2)

% 3) e1, stable between R=R_tra1 and R= R_tra2;
plot3(R_zoom(1:5),e1_zoom(1,1:5),...
    mult*e1_zoom(2,1:5),'--k','LineWidth',1.5)
plot3(R_zoom(5:end),e1_zoom(1,5:end),...
    mult*e1_zoom(2,5:end),'-k','LineWidth',1.5)


% 4) e2, stable between  R=0 and R=R_tra1
plot3(R_zoom(1:5),e2_zoom(1,1:5),...
    mult*e2_zoom(2,1:5),'-k','LineWidth',1.5)
plot3(R_zoom(5:end),e2_zoom(1,5:end),...
    mult*e2_zoom(2,5:end),'--k','LineWidth',1.5)


xlim([0 0.015])
ylim([0 0.06])
zlim([0 mult*0.0001])
xticks([0 0.0075 0.015])
yticks([0 0.03 0.06])
zticks([0 0.05 0.1])
xticklabels([])
yticklabels([])
zticklabels([])
set(gca,'FontSize',15)
view([-70,40])
grid on

hold off
%%

% 8) the 2 parameter diagram

GH_ind = find(abs(H2par(:,1)-GH(1))<0.005);
figure
hold on
plot(...
    R_tra1_2par,del1,...
    R_tra2_2par,del1,...
    'LineWidth',2,'Color',[0.4 0.83 0.1])
plot(...
    H2par(1:GH_ind,1),H2par(1:GH_ind,2),'r--',...
    H2par(GH_ind:end,1),H2par(GH_ind:end,2),'r-',...
    h2par(1:Rhind,1),h2par(1:Rhind,2),'--k',...
    h2par(Rhind:end,1),h2par(Rhind:end,2),'-k',...
    'LineWidth',2)
plot(...
     Fl2par(:,1),Fl2par(:,2),'Color',[0.07,0.6,1],...
    'LineWidth',2)

plot(TT(1),TT(2),'.k','MarkerSize',24)
plot(GH(1),GH(2),'.k','MarkerSize',24)
plot(Rh(1),Rh(2),'.k','MarkerSize',24)
plot(linspace(0,8,10),2.4*ones(1,10),'--','Color',[.7 .7 .7]);

xlabel('$R$')
ylabel('$\delta$','Rotation',0)
axis([0 3 0 2.7])
box on
set(gca,'FontSize',14)

% bifurcation zoomed in

axes('Position',[0.55 0.2 0.3 0.4])
hold on

plot(...
    R_tra1_2par,del1,...
    R_tra2_2par,del1,...
    'LineWidth',1.5,'Color',[0.4 0.83 0.1])

plot(...
    H2par(1:GH_ind,1),H2par(1:GH_ind,2),'r--',...
    H2par(GH_ind:end,1),H2par(GH_ind:end,2),'r-',...
    h2par(1:Rhind,1),h2par(1:Rhind,2),'--k',...
    h2par(Rhind:end,1),h2par(Rhind:end,2),'-k',...
    'LineWidth',1.5)

plot(...
     Fl2par(:,1),Fl2par(:,2),'Color',[0.07,0.6,1],...
    'LineWidth',1.5)

plot(GH(1),GH(2),'ok',...
    Rh(1),Rh(2),'ok',...
    'LineWidth',3.5,'MarkerSize',3)

axis([0.15 1.35 .3 1.54])
box on
set(gca,'XTick',[0.15 1.35], 'YTick', [.3 1.54])

annotation('textbox',...
    [0.639 0.825 0.0761 0.077],...
    'String','$H$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','on');
annotation('textbox',...
    [0.4767 0.6440 0.0648 0.0771],...
    'String',{'$h$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','on');
annotation('textbox',...
    [0.15 0.725 0.0784 0.077],...
    'String','$T_1$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14);
annotation('textbox',...
    [0.3482 0.749 0.0784 0.0771],...
    'String','$T_2$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14);

annotation('textbox',...
    [0.63 0.237 0.094 0.077],...
    'String','$GH$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','on');

set(gca,'FontSize',14)
%% phase portraits

% Rphas1  = 2.47;
% Rphas2  = 3;
% odefun1  = @(t,var)RModefun(var,Rphas1);
% odefun2  = @(t,var)RModefun(var,Rphas2);
% perinit1 = [5.1016 0.0074];
% perinit2 = [11.9229 0.01533];
% 
% T1 = 8.5;
% T2 = 15;
% 
% per_phas1    = ode45(odefun1,[0,T1],perinit1);
% per_phas2    = ode45(odefun2,[0,T2],perinit2);
% varper_phas1 = deval(per_phas1,linspace(0,T1,250));
% varper_phas2 = deval(per_phas2,linspace(0,T2,250));
% 
% e0_phas1    = [0,0];
% e1_phas1    = [Rphas1/C, 0];
% e2_phas1    = [mu,0];
% e3_phas1(1) = (delta*beta)/((gamma*alpha) - delta);
% e3_phas1(2) = (Rphas1-(C*e3_phas1(1)))*(beta+e3_phas1(1))*(e3_phas1(1)-mu)...
%     /(alpha*(nu + e3_phas1(1)));
% 
% % the stable manifold of e1
% maninit1 = e1_phas1 + [-0.01,0];
% maninit2 = e1_phas1 + [0.01,0];
% maninit3 = e2_phas1 + [0 1e-8];
% 
% [~,varman_phas1_1]   = ode45(odefun1,[15,0],maninit1,opts);
% [~,varman_phas1_2]   = ode45(odefun1,[15,0],maninit2,opts);
% [~,varman_phas1_3]   = ode45(odefun1,[15,0],maninit3,opts);
% 
% figure;
% hold on
% box on
% 
% plot(...
%     e0_phas1(1),e0_phas1(2),'o',...
%     'MarkerSize',3,'LineWidth',3.5,'Color','k')
% plot(...
%     e1_phas1(1),e1_phas1(2),'+',...
%     e2_phas1(1),e2_phas1(2),'+',...
%     'MarkerSize',7,'LineWidth',1.5,'Color','k')  
% plot(...
%     e3_phas1(1),e3_phas1(2),'o',...
%     'MarkerSize',5,'LineWidth',2,'Color','k')
% 
% plot(varper_phas1(1,:),varper_phas1(2,:),'-k','LineWidth',2)
% 
% plot(varper_phas2(1,:),varper_phas2(2,:),'Color',[.07 .6 1],'LineWidth',2)
% 
% plot(...
%     varman_phas1_1(:,1),varman_phas1_1(:,2),'--k',...
%     varman_phas1_2(:,1),varman_phas1_2(:,2),'--k',...
%     varman_phas1_3(:,1),varman_phas1_3(:,2),'--k')
% 
% axis([-0.5 16 -.001 0.03])
% set(gca,'FontSize',14)
% xlabel('$R$')
% ylabel('$P$','Rotation',0)
% 
%% functions

function [dvar] = RModefun(var,RR)
%parameters
global C gamma beta alpha mu nu

delta = 2.2;

N = var(1);
P = var(2);

dN = RR * N *(1-((C/RR)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = gamma * ((alpha*P*N)/(beta + N)) - (delta*P);

dvar = [dN;dP];
end
