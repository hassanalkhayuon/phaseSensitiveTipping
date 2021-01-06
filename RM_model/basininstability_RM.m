% to fined the limit cycle gamma(t,R) for each R_H<R<R_h by a BVP

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');


% parameters


R_star     = 2.5;
del_star   = 2.2;    %/yr.            Predator death rate in absent of prey
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter

%% step 1 defining the initial conitions at (R_star,del_star)

par = [R_star,del_star];
initcond = [2.8    0.007];
Tper     = 12;
T        = 1000;
ivpfun   = @(t,var)RModefun(var,par);
perfun   = @(t,var,Tper)Tper*RModefun(var,par);

[~,tempvar]  = ode45(ivpfun,[0 T],initcond);
initper      = [3.300060384770406,0.022423830878037];
ivpsol       = ode45(@(t,var)perfun(t,var,Tper),[0 1],initper);


% solinit
ss = linspace(0,1,200);
tempinit = @(s)deval(ivpsol,s);
solinit=bvpinit(ss,tempinit,Tper);

% The boundary condtions for the periodic orbit
BC=@(WL,WR,Tper)(...
   [WL(1)-WR(1);...
    WL(2)-WR(2);...
    WL(2)-initper(2);... %point phase condition
    ]);

BVPper_start= bvp5c(perfun,BC,solinit);


N_initconds = 200; % number of initial conditions on the limit cycle

% Tend1   = 0.9539875;
% Tend2   = 0.0480805;
% 
% ss1 = 1:-0.01:Tend1;
% ss2 = 0:0.01:Tend2;
% ss = [flip(ss1),(ss2)];

ss = linspace(0,1,N_initconds);
% load('perinitstart_RM.mat')
% initconds = deval(ivpsol,ss);
initconds = deval(BVPper_start,ss);
N_initconds = length(initconds);
figure
plot(initconds(1,:),1000*initconds(2,:),'-')
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%% Step 2 a grid in the oscillatory coexistence region

load('RM2par_new.mat')

ind1H    = 935;
ind2H    = 1033;
ind1h    = 110;
ind2h    = 453;


D        = 4;
L1       = floor((ind2H - ind1H)/D);
L2       = floor((ind2h - ind1h)/D);

RH_bif   = NaN(1,L1);
delH_bif = NaN(1,L1);

Rh_bif   = NaN(1,L2);
delh_bif = NaN(1,L2);



n=1;
for ind_bif =ind1H:ind2H
    if rem(ind_bif,D)==0
        RH_bif(n)   = H2par(ind_bif,1);
        delH_bif(n) = H2par(ind_bif,2);
        n=n+1;
    end
end

n=1;
for ind_bif =ind1h:ind2h
    if rem(ind_bif,D)==0
        Rh_bif(n)   = h2par(ind_bif,1);
        delh_bif(n) = h2par(ind_bif,2);
        n=n+1;
    end
end


delH = @(RR)interp1(RH_bif,delH_bif,RR);
delh = @(RR)interp1(Rh_bif,delh_bif,RR);

Rmin = 0.35;
Rmax = 2.47;

delmax = 2.7;
delmin = 0.45;

T = 100; % integration time

grid_R   = 500;
grid_del = grid_R;


R_scan   = linspace(Rmin,Rmax,grid_R);
del_scan = linspace(delmin,delmax,grid_del);
tip = zeros(grid_del,grid_R);
tic

for ind_R = 1:grid_R
    R = R_scan(ind_R);
    if and(R>Rmin,R<Rmax)
        del1 = delh(R);
        del2 = delH(R);
        for ind_del = 1:grid_del
            del = del_scan(ind_del);
            if and(del1<del,del<del2)
                %                 tip(ind_del,ind_R)=1;
                par = [R,del];
                odefun = @(t,var)RModefun(var,par);
                for ind_initcond = 1: N_initconds
                    [~,var]=ode23tb(odefun,[0 T],initconds(:,ind_initcond));
                    if norm(var(end,:))<1e-3
                        tip(ind_del,ind_R) = tip(ind_del,ind_R) + 1;
                    end
                end
            end
             disp([ind_R,ind_del]);
        end
    end
end

%%
initpoint = initconds;
R_thre = linspace(1.68,2.356,grid_R);
del_thre = NaN(size(R_thre));
del_init = 2.28;
options = optimset('TolX',1e-6);
for ind_thre = 1:grid_R
    R = R_thre(ind_thre);
    GG = @(del)tipthre(R,del,initpoint,T);
    del_thre(ind_thre) = fzero(GG,del_init);
    del_init = del_thre(ind_thre);
    disp(ind_thre);
end

%%
grayscal = 1:-0.01:0;
map = [grayscal',grayscal',grayscal'];
% figure;hold on
PCOLOR = pcolor(R_scan,del_scan,tip);
PCOLOR.LineStyle = 'none';
PCOLOR.FaceAlpha = 0.2;
colormap(map)
caxis([0 1]);
% plot(RH_bif,delH_bif,'r-')
% plot(Rh_bif,delh_bif,'r-')
% plot(R_star,del_star,'k.','MarkerSize',20)
plot(R_thre,del_thre,'-','LineWidth',1.5,'Color',[.65 .65 .65])


%% phase portraits

R   = 2.313;
del = 2.123;
GG = @(del)tipthre(R,del,initpoint,T);
GG(del)
par = [R,del];
odefun  = @(t,var)RModefun(var,par);

e0    = [0,0];
e1   = [R/C, 0];
e2    = [mu,0];
e3(1) = (del*beta)/((gamma*alpha) - del);
e3(2) = (R-(C*e3(1)))*(beta+e3(1))*(e3(1)-mu)...
    /(alpha*(nu + e3(1)));

G   = @(var)RModefun(var,par);
JJ = MyJacobian(G,e1);

[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert1 = eigvic(:,2);
else
    pert1 = eigvic(:,1);
end
pert1 = 1e-6*pert1';

G   = @(var)RModefun(var,par);
JJ = MyJacobian(G,e2);

[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert2 = eigvic(:,2);
else
    pert2 = eigvic(:,1);
end
pert2 = 1e-6*pert2';

% % the stable manifold of e1
maninit11 = e1 + pert1;
maninit12 = e1 - pert1;
maninit2  = e2 + pert2;

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~,varman1]   = ode45(odefun,[50,0],maninit11,opts);
[~,varman2]   = ode45(odefun,[50,0],maninit12,opts);
[~,varman3]   = ode45(odefun,[50,0],maninit2,opts);

 figure;
hold on
box on

plot(...
    e0(1),e0(2),'o',...
    'MarkerSize',3,'LineWidth',3.5,'Color','k')
plot(...
    e1(1),e1(2),'+',...
    'MarkerSize',7,'LineWidth',1.5,'Color','k')
plot(...
    e3(1),e3(2),'o',...
    'MarkerSize',5,'LineWidth',2,'Color','k')

plot(...
    varman1(:,1),varman1(:,2),'--k',...
    varman2(:,1),varman2(:,2),'--k',...
    varman3(:,1),varman3(:,2),'--k')

plot(initconds(1,:),initconds(2,:),'-b')

axis([-0.5 20 -.001 0.04])
set(gca,'FontSize',14)
xlabel('$R$')
ylabel('$P$','Rotation',0)

%% functions

function [crit] = tipthre(R,del,initpoints,T)
par = [R,del];
parfor ind = 1:length(initpoints)
    initpoint = initpoints(:,ind);
    odefun=@(t,var)RModefun(var,par);
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [~,varshot] = ode45(odefun,[0,T],initpoint,opts);
    ccc(ind) = norm(varshot(end,:))
end
if min(ccc) < 1e-4
    crit = 1;
else
    crit  = -1;
end
end

%the ode fuction of RM model
