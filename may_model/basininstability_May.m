% to analyse basin-instability of the May Model.

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');



% parameters:

R_star  = 3.3;       % /yr.              Prey intrinsic growth rate
q_star  = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

 T        = 100;
 
 opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
%% step 1 defining the initial conditios at (R_start,q_start)

par = [R_star,q_star];
initcond = [5.3   0.02183];
Tper     = 9.5;

ivpfun   = @(t,var)Mayodefun(var,par);
perfun   = @(t,var,Tper)Tper*Mayodefun(var,par);

[~,tempvar]  = ode45(ivpfun,[0 T],initcond);
initper  = tempvar(end,:); %the top point is [6.9261, 0.0339];
ivpsol = ode45(@(t,var)perfun(t,var,Tper),[0 1],initper);


% The boundary condtions for the periodic orbit
ss = linspace(0,1,1000);
tempinit = @(s)deval(ivpsol,s);
solinit=bvpinit(ss,tempinit,Tper);

BC=@(WL,WR,Tper)(...
    [WL(1)-WR(1);...
    WL(2)-WR(2);...
    WL(2)-initper(2);... %point phase condition
    ]);

BVPper_start= bvp5c(perfun,BC,solinit);

N_initconds = 50; % number of initial conditions on the limit cycle
ss = linspace(0,1,N_initconds);
% load('perinitstart_may.mat')
initconds = deval(BVPper_start,ss);
figure(3); hold on
plot(initconds(1,:),initconds(2,:))
opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
%% step 2 a grid in the oscillatory coexistence region

load('may2par_xppaut.mat')

D = 20;
ind1H1 = 700;
ind2H1 = 1645;

ind1H2 = 1645;
ind2H2 = 2060;

ind1H3 = 2000;
ind2H3 = 3300;

ind1h  = 680;
ind2h  = 1200;

L1 = floor((ind2H1 - ind1H1)/D);
L2 = floor((ind2H2 - ind1H2)/D);
L3 = floor((ind2H3 - ind1H3)/D);
L4 = floor((ind2h - ind1h)/D);

RH1_bif=NaN(1,L1);
qH1_bif=NaN(1,L1);

RH2_bif=NaN(1,L2);
qH2_bif=NaN(1,L2);

RH3_bif=NaN(1,L3);
qH3_bif=NaN(1,L3);

Rh_bif=NaN(1,L4);
qh_bif=NaN(1,L4);

n=1; % H1
for ind_bif =ind1H1:ind2H1
    if rem(ind_bif,D)==0
        RH1_bif(n)=H2par(ind_bif,1);
        qH1_bif(n)=H2par(ind_bif,2);
        n=n+1;
    end
end

n=1; % H2
for ind_bif =ind1H2:ind2H2
    if rem(ind_bif,D)==0
        RH2_bif(n)=H2par(ind_bif,1);
        qH2_bif(n)=H2par(ind_bif,2);
        n=n+1;
    end
end

n=1; %H3
for ind_bif =ind1H3:ind2H3
    if rem(ind_bif,D)==0
        RH3_bif(n)=H2par(ind_bif,1);
        qH3_bif(n)=H2par(ind_bif,2);
        n=n+1;
    end
end

n=1; %h
for ind_bif =ind1h:ind2h
    if rem(ind_bif,D)==0
        Rh_bif(n)=h2par(ind_bif,1);
        qh_bif(n)=h2par(ind_bif,2);
        n=n+1;
    end
end

qH1 = @(RR)interp1(RH1_bif,qH1_bif,RR);
qH2 = @(RR)interp1(RH2_bif,qH2_bif,RR);
qH3 = @(RR)interp1(RH3_bif,qH3_bif,RR);
qh  = @(RR)interp1(Rh_bif,qh_bif,RR);

R1 = 1.549;
R2 = 1.69;
R3 = 4;

T = 100; %integration time

grid_q = 200;
grid_R = 2*grid_q;
R_scan = linspace(1.55,3.7,grid_R);
q_scan = linspace(160,320,grid_q);
tip = zeros(grid_q,grid_R);
tic
for ind_R = 1:grid_R
    R = R_scan(ind_R);
    if and(R>=R1,R<=R2)
        q1 = qH1(R);
        q2 = qH2(R);
        parfor ind_q = 1:grid_q
            q = q_scan(ind_q);
            if and(q1<q,q<q2)
%                 tip(ind_q,ind_R)=1;
                par = [R,q];
                odefun = @(t,var)Mayodefun(var,par);
                for ind_initcond = 1: N_initconds
                    [~,var]=ode45(odefun,[0 T],initconds(:,ind_initcond));
                    e0 = [0,eps/q];
                    if norm(var(end,:)-e0)<1e-3
                        tip(ind_q,ind_R) = tip(ind_q,ind_R) + 1;
                    end
                end
            end
        end
    elseif and(R>R2,R<=R3)
        q1 = qh(R);
        q2 = qH3(R);
        parfor ind_q = 1:grid_q
            q = q_scan(ind_q);
            if and(q1<q,q<q2)
%                 tip(ind_q,ind_R)=1;
                par = [R,q];
                odefun = @(t,var)Mayodefun(var,par);
                for ind_initcond = 1: N_initconds
                    [~,var]=ode45(odefun,[0 T],initconds(:,ind_initcond));
                    e0 = [0,eps/q];
                    if norm(var(end,:)-e0)<1e-3
                        tip(ind_q,ind_R) = tip(ind_q,ind_R) + 1;
                    end
                end
            end
        end
    end
    disp(ind_R);
end

%% the threshold
sper_max = 0.785578557855786;
initpoint = deval(BVPper_start,sper_max);
R_thre = linspace(1.57,3.3,200);
q_thre = NaN(size(R_thre));
q_init = 260;
options = optimset('TolX',1e-6);
for ind_thre = 1:200
    R = R_thre(ind_thre);
    GG = @(qq)tipthre(R,qq,initpoint,T);
    q_thre(ind_thre) = fzero(GG,q_init);
    q_init = q_thre(ind_thre);
    disp(ind_thre);
end
%%
grayscal = 1:-0.01:0;
map = [grayscal',grayscal',grayscal'];
figure;hold on
PCOLOR = pcolor(R_scan,q_scan,tip);
PCOLOR.FaceAlpha = 0.5;
PCOLOR.LineStyle = 'none';
colormap(map)
caxis([0 N_initconds]);
plot(RH2_bif,qH2_bif,'r-')
plot(RH1_bif,qH1_bif,'r-')
plot(RH3_bif,qH3_bif,'r-')
plot(Rh_bif,qh_bif,'r-')
plot(R_star,q_star,'m.','MarkerSize',20)
plot(R_thre,q_thre,'-','LineWidth',1.5,'Color',[.8 .8 .8])

% 4) (R,q) parameter plane
load('may2par_xppaut.mat')

BTind = 2751;
GHind = 973;
Rhind = 963;

BT = [1.77476, 137.231];
GH = [1.64885, 209.740];
Rh  = h2par(Rhind,:);

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

plot(...
    BT(1),BT(2),'k.',...
    GH(1),GH(2),'k.',...
    Rh(1),Rh(2),'k.',...
    'MarkerSize',20)
% plot((C*mu)*ones(1,100),linspace(0,450),'Color',[0.4 0.83 0.1],'LineWidth',1.5)


axis([1 3.7 100 320])
xlabel('$R$')
ylabel('$q$','Rotation',0)
set(gca,'FontSize',16)
box on


%% phase portraits

R_crit = 3.3;
q_crit = 205;
par1   = [R_crit,q_crit];
odefun1 = @(t,var)Mayodefun(var,par1);
% equilibria

[N_eq] = May_eq(R_crit,q_crit);

e0    = [0,eps/q_crit];
e1    = [par1(1)/C, 0];
e2    = [mu,0];

e3(1) =  N_eq(1);
e3(2) = (N_eq(1) + eps)/q_crit;

e4(1) =  N_eq(2);
e4(2) = (N_eq(2) + eps)/q_crit;

e5(1) =  N_eq(3);
e5(2) = (N_eq(3) + eps)/q_crit;

G   = @(var)Mayodefun(var,par1);
JJ = MyJacobian(G,e5);

[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-4*pert';

maninit1 = e5 + pert;
maninit2 = e5 - pert;

[~,varman1]   = ode45(odefun1,[11,0],maninit1,opts);
[~,varman2]   = ode45(odefun1,[11,0],maninit2,opts);
figure;
hold on
plot(initconds(1,:),initconds(2,:),'-b')
plot(...
    e0(1),e0(2),'ok',...
    'MarkerSize',4,'LineWidth',3.5)
plot(...
    e1(1),e1(2),'+k',...
    e5(1),e5(2),'+k',...
    'MarkerSize',6,'LineWidth',1)
plot(...
    e2(1),e2(2),'ok',...
    e3(1),e3(2),'ok',...
    e4(1),e4(2),'ok',...
    'MarkerSize',6,'LineWidth',1)

plot(...
    varman1(:,1),varman1(:,2),'--k',...
    varman2(:,1),varman2(:,2),'--k')

axis([-0.5 16 -.001 0.05])
set(gca,'FontSize',14)
xlabel('$R$')
ylabel('$P$','Rotation',0)


%% Functions

% the flip function
function [crit] = tipthre(R,q,initpoint,T)
e0 = [0,eps/q];
par = [R,q];
odefun=@(t,var)Mayodefun(var,par);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~,varshot] = ode45(odefun,[0,T],initpoint,opts);
if norm(varshot(end,:)-e0) < 1e-3
    crit = 1;
else
    crit  = -1;
end
end

% the ode function 
function [dvar] = Mayodefun(var,par)

%parameters

C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.


N = var(1);
P = var(2);

RR = par(1);
q  = par(2);

dN = RR * N *(1-((C/RR)*N))*((N - mu)/(nu + N)) - ((alpha*P*N)/(beta + N));
dP = s*P*( 1-((q*P)/(N+eps)) );

dvar = [dN;dP];

end

