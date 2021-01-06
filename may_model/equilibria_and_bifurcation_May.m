% the equlibria and stability of May model.

% for detailes about the model and parameters go to:
% https://www.overleaf.com/project/5db3114c13676a00018d4902

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
clear
clc

% parameters:
global  C alpha beta s mu nu eps

R       = 0:0.001:4.5;     %/yr.          Prey intrinsic growth rate
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    %/yr.               Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr                  Rate of predator population.
mu      = 0.03;      % NA                 Allee parameter.
nu      = 0.003;      % NA                Allee parameter.
eps     = 0.031;       % NA                Artificial parameter.

R_zoom  = [linspace(0,C*mu,5),linspace(C*mu+0.001,0.015,5)];

% poly = May_eq(4)
% roots(poly)
%%

e0 = NaN(2,length(R));
e1 = NaN(2,length(R));
e2 = NaN(2,length(R));
e3 = NaN(2,length(R));
e4 = NaN(2,length(R));
e5 = NaN(2,length(R));

e0_zoom = NaN(2,length(R_zoom));
e1_zoom = NaN(2,length(R_zoom));
e2_zoom = NaN(2,length(R_zoom));

for nn = 1:length(R)
    RR = R(nn);
    
    e0(1,nn) = 0;
    e0(2,nn) = eps/q;
    
    e1(1,nn) = RR/C;
    e1(2,nn) = 0;
    
    e2(1,nn) = mu;
    e2(2,nn) = 0;
    
    [N_eq] = May_eq(RR,q);
    
    e3(1,nn) =  N_eq(1);
    e3(2,nn) = (N_eq(1) + eps)/q;
    
    e4(1,nn) =  N_eq(2);
    e4(2,nn) = (N_eq(2) + eps)/q;
    
    e5(1,nn) =  N_eq(3);
    e5(2,nn) = (N_eq(3) + eps)/q;
    
    
    % stability
    par = [RR,q]; 
    G   = @(var)Mayodefun(var,par);
    
    %e0
    JJ = MyJacobian(G,e0(:,nn));
    if isfinite(sum(sum(JJ)))
        stability0(:,nn)=eig(JJ);
    else
        stability0(:,nn)=NaN(2,1);
    end
    
    %e1
    JJ = MyJacobian(G,e1(:,nn));
    if isfinite(sum(sum(JJ)))
        stability1(:,nn)=eig(JJ);
    else
        stability1(:,nn)=NaN(2,1);
    end
    
    %e2
    JJ = MyJacobian(G,e2(:,nn));
    if isfinite(sum(sum(JJ)))
        stability2(:,nn)=eig(JJ);
    else
        stability2(:,nn)=NaN(2,1);
    end
    
    %e3
    JJ = MyJacobian(G,e3(:,nn));
    if isfinite(sum(sum(JJ)))
        stability3(:,nn)=eig(JJ);
    else
        stability3(:,nn)=NaN(2,1);
    end
    
    %e4
    JJ = MyJacobian(G,e4(:,nn));
    if isfinite(sum(sum(JJ)))
        stability4(:,nn)=eig(JJ);
    else
        stability4(:,nn)=NaN(2,1);
    end
    
    %e5
    JJ = MyJacobian(G,e5(:,nn));
    if isfinite(sum(sum(JJ)))
        stability5(:,nn)=eig(JJ);
    else
        stability5(:,nn)=NaN(2,1);
    end
end

for mm =1:length(R_zoom)
    RR = R_zoom(mm);
    
    e0_zoom(1,mm) = 0;
    e0_zoom(2,mm) = eps/q;
    
    e1_zoom(1,mm) = RR/C;
    e1_zoom(2,mm) = 0;
    
    e2_zoom(1,mm) = mu;
    e2_zoom(2,mm) = 0;
end
%% bifurcation points (xppaut)
R_HB1 = 1.662;
R_HB2 = 3.81;
R_SN  = 1.205;
R_Fl  = 1.661;
HB_ind1 = find(abs(R-R_HB1)<0.003);
HB_ind2 = find(abs(R-R_HB2)<0.0001);
SN_ind  = find(abs(R-R_SN)<0.001);


%% plots
set(0,'defaulttextInterpreter','latex')

mult = 1000;

figure
hold on

% 1) bifurcation points
plot3(...
    R(SN_ind),e4(1,SN_ind),mult*e4(2,SN_ind),'k.',... % SN point
    R(HB_ind1),e4(1,HB_ind1),mult*e4(2,HB_ind1),'k.',...% first HB point
    R(HB_ind2),e3(1,HB_ind2),mult*e3(2,HB_ind2),'k.',...% second HB point
    C*mu,mu,0,'k.',...                               % transcritical point
    'MarkerSize',20)

plot3(...
    R,e0(1,:),mult*e0(2,:),'k-',...
    R,e1(1,:),mult*e1(2,:),'k--',...
    R,e2(1,:),mult*e2(2,:),'k--',...
    R(1:HB_ind1),e4(1,1:HB_ind1),mult*e4(2,1:HB_ind1),'k-',...
    R(HB_ind1:end),e4(1,HB_ind1:end),mult*e4(2,HB_ind1:end),'k--',...
    R(1:HB_ind2),e3(1,1:HB_ind2),mult*e3(2,1:HB_ind2),'k--',... 
    R(HB_ind2:end),e3(1,HB_ind2:end),mult*e3(2,HB_ind2:end),'k-',...
    R,e5(1,:),e5(2,:),'k--',...
    'LineWidth',1.5)

box on
grid on
view([-70 40])
load('data\may_per.mat')
for kk = 1:length(Rper)
    ss=linspace(0,1,500);
    varper=deval(per(kk),ss);
    plot3(Rper(kk)*ones(size(ss)),varper(1,:),mult*varper(2,:),...
        'k','LineWidth',1)
end
for mm = 1:3
    ss=linspace(0,1,500);
    varper_un=deval(per_un(mm),ss);
    plot3(Rper_un(mm)*ones(size(ss)),varper_un(1,:),mult*varper_un(2,:),...
        '--k','LineWidth',1) 
end
axis([0 4.5 0 25 0 0.04*mult])

xlabel('$R$')
ylabel('$N$')
zlabel('$P \times 10^3$')
set(gca,'FontSize',14)


%% 3) Zoomed in
axes('Position',[0.252 0.65 0.25 0.25])
box on
hold on
grid on
plot3(...
    R,e5(1,:),mult*e5(2,:),'k--',...
    R(1:HB_ind1),e4(1,1:HB_ind1),mult*e4(2,1:HB_ind1),'k-',...
    R(HB_ind1:end),e4(1,HB_ind1:end),mult*e4(2,HB_ind1:end),'k--',...
    'LineWidth',1.5)
plot3(...
    R(SN_ind),e4(1,SN_ind),mult*e4(2,SN_ind),'k.',... % SN po
    'MarkerSize',24)

axis([1.2 1.7 0 .4 0 1.5])
xticks([1.2 1.45 1.7])
yticks([0 .2 .4])
zticks([0 .75 1.5])
xticklabels([])
yticklabels([])
zticklabels([])
view([-70 40])

axes('Position',[0.15 0.255 0.25 0.25])
box on
hold on

plot3(...
    R_zoom,e0_zoom(1,:),mult*e0_zoom(2,:),'k-',...
    R_zoom(1:5),e1_zoom(1,1:5),mult*e1_zoom(2,1:5),'k.',...
    R_zoom(5:end),e1_zoom(1,5:end),mult*e1_zoom(2,5:end),'k--',...
    R_zoom(1:5),e2_zoom(1,1:5),mult*e2_zoom(2,1:5),'k--',...
    R_zoom(5:end),e2_zoom(1,5:end),mult*e2_zoom(2,5:end),'k:',...
    'LineWidth',1.5)
plot3(C*mu,mu,0,'k.','MarkerSize',24)

axis([0 0.015 0 0.07 0 0.2])
xticks([0 0.0075 0.015])
yticks([0 0.035 0.07])
zticks([0 0.1 0.2])
xticklabels([])
yticklabels([])
zticklabels([])
box on
grid on
view([-70 40])

%%
% 4) (R,q) parameter plane
load('data\may2par_xppaut.mat')

BTind = 2751;
GHind = 973;
Rhind = 963;

BT = [1.77476, 137.231];
GH = [1.64885, 209.740];
Rh  = h2par(Rhind,:);




figure;
hold on
plot(...
    F2par(1:BTind,1),F2par(1:BTind,2),'-',...
    F2par(BTind:end,1),F2par(BTind:end,2),'--',...
    'LineWidth',1.5,'Color','b');
plot(...
    H2par(1:GHind,1),H2par(1:GHind,2),'--r',...
    H2par(GHind:end,1),H2par(GHind:end,2),'-r',...
    h2par(1:Rhind,1),h2par(1:Rhind,2),'--k',...
    h2par(Rhind:end,1),h2par(Rhind:end,2),'-k',...
    'LineWidth',1.5);

plot(Fl22par(:,1),Fl22par(:,2),'-','Color',[.07 .6 1],'LineWidth',1.5);

plot(linspace(0,10),q*ones(1,100),'--','Color',[.8 .8 .8])
plot(...
    BT(1),BT(2),'k.',...
    GH(1),GH(2),'k.',...
    Rh(1),Rh(2),'k.',...
    'MarkerSize',20)
plot((C*mu)*ones(1,100),linspace(0,450),'Color',[0.4 0.83 0.1],'LineWidth',1.5)


axis([0 4 50 350])
xlabel('$R$')
ylabel('$q$','Rotation',0)
set(gca,'FontSize',16)
box on

%% scan 

R_scan   =0.00000001;
q_scan   = 128.7;
par_scan = [R_scan,q_scan];
odefun   = @(t,var)Mayodefun(var,par_scan);
opts     = odeset('RelTol',1e-4,'AbsTol',1e-8);

% equilibria 
e0_scan(1) = 0;
e0_scan(2) = eps/q_scan;

e1_scan(1) = R_scan/C;
e1_scan(2) = 0;

e2_scan(1) = mu;
e2_scan(2) = 0;

[N_eq] = May_eq(R_scan,q_scan);

e3_scan(1) =  N_eq(1);
e3_scan(2) = (N_eq(1) + eps)/q_scan;

e4_scan(1) =  N_eq(2);
e4_scan(2) = (N_eq(2) + eps)/q_scan;

e5_scan(1) =  N_eq(3);
e5_scan(2) = (N_eq(3) + eps)/q_scan;

FF   = @(var)Mayodefun(var,par_scan);
JJ   = MyJacobian(FF,e1_scan);

[eigvic,eigval]=eig(JJ);
if eigval(2,2)>0
    pert1 = eigvic(:,2);
else
    pert1 = eigvic(:,1);
end
pert1 = 1e-4*pert1';

% maninit1 = e5 + pert;
% maninit2 = e5 - pert;
maninit3 = e1_scan + pert1;
maninit4 = e1_scan - pert1;

% [~,varman1]   = ode45(odefun,[20,0],maninit1,opts);
% [~,varman2]   = ode45(odefun,[20,0],maninit2,opts);
[~,varman3]   = ode45(odefun,[0,100],maninit3,opts);
[~,varman4]   = ode45(odefun,[0,100],maninit4,opts);

scangrid = 10;
N_end    = 10;
P_end    = 0.04;
N_scan   = linspace(0.01,N_end,scangrid);
P_scan   = linspace(0.00001,P_end,scangrid);

T        = 100;

figure;
hold on
plot(...
    e2_scan(1),e2_scan(2),'+k',...
    'MarkerSize',8,'LineWidth',1)
plot(...
    e1_scan(1),e1_scan(2),'ok',...
    'MarkerSize',3,'LineWidth',1)
plot(...
    varman3(:,1),1000*varman3(:,2),'r--',...
    varman4(:,1),varman4(:,2),'r--')

axis([0 N_end 0 1000*P_end]);
box on

for ind_scanN     = 1:scangrid
    for ind_scanP = 1:scangrid
        initcond  = [N_scan(ind_scanN),P_scan(ind_scanP)];
        [~,var]   = ode45(odefun,[0 T],initcond,opts);
        
        plot(var(:,1),1000*var(:,2),'k')
        drawnow;
    end
end
%%
function [dvar] = Mayodefun(var,par)
%parameters
global  C alpha beta s mu nu eps
N = var(1);
P = var(2);

RR = par(1);
q  = par(2);

dN = RR * N *(1-((C/RR)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = s*P*( 1-((q*P)/(N+eps)) );

dvar = [dN;dP];
end