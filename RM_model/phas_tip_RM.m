% phas_tip_RM
% To find the tipping phases of the RM model.
% the starting point is at R = 2.47, del = 2.2. 
% del is going to be fixes through out, and R is going to be decreased up
% to R_H = 1.67.

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

R_star     = 2.5;
R_end      = 1.6;
del        = 2.2;    %/yr.            Predator death rate in absent of prey
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter

%% initiation 

N_initconds = 5000; % number of initial conditions on the limit cycle

% Tend1   = 0.9539875;
% Tend2   = 0.0480805;

Tend1   = 0.97;
Tend2   = 0.03; % for the red curve 


ss1 = 1:-0.0001:Tend1;
ss2 = 0:0.0001:Tend2;
ss = [flip(ss1),(ss2)];
ss  = linspace(0,1,N_initconds);

% load('perinitstart_RM.mat')     % R=2.47
load('perinitstart_RM1.mat')      % R = 2.5
varper = deval(BVPper_start,ss);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

e3_star(1) = (del*beta)/((gamma*alpha) - del);
e3_star(2) = (R_star-(C*e3_star(1)))*(beta+e3_star(1))*(e3_star(1)-mu)...
    /(alpha*(nu + e3_star(1)));

[thper,rper]  = cart2pol(varper(1,:)-e3_star(1),(1000*(varper(2,:)-e3_star(2))));
AA = sortrows([thper',rper']);
thper = AA(:,1);
rper  = AA(:,2);
Rscan = linspace(R_end,1.95,10);

Tintg  = 100;

figure(1); hold on
plot(rper.*cos(thper)+e3_star(1),1000*((rper.*sin(thper)/1000)+e3_star(2)),'r-','LineWidth',2)
% plot(e3_star(1),e3_star(2),'.k','MarkerSize',15)
%% scan

for indR = 1:length(Rscan)
    R = Rscan(indR);
    TIPind = 1;
    parfor indth = 1:length(thper)
        N = rper(indth)*cos(thper(indth))+ e3_star(1);
        P = (rper(indth)*sin(thper(indth))/1000)+ e3_star(2);
        init    = [N;P];
        par     = [R,del];
        odefun  = @(t,var)RModefun(var,par);
        [t,var] = ode45(odefun,[0 Tintg],init,opts);
        if norm(var(end,:))<=1e-4
            tip(indth,indR) = 1;
            TIPind = TIPind+1;
        else
            tip(indth,indR) = 0;
        end
    end
    tip(:,indR)= tip(:,indR) * TIPind;
    disp([indR,length(Rscan)])
end





%% threshold
Ttop    = 1;
Tend1   = 0.9539875;
Tend2   = 0.0480805;
Tinit1  = 1:-0.002:Tend1;
Tinit2  = 0:0.001:Tend2;
Tintg   = 100;
Rinit   = R_end;

Rcrit1  = NaN(size(Tinit1));
thcrit1 = NaN(size(Tinit1));
rcrit1  = NaN(size(Tinit1));

for ind_T = 1:length(Tinit1)
    initpoint        = deval(BVPper_start,Tinit1(ind_T));
    GG               = @(RR)tipthre(RR,del,initpoint,Tintg);
    Rinit            = fzero(GG,Rinit);
    Rcrit1(ind_T)    = Rinit;
    Rinit            = 1.6;
    x                = initpoint(1)-e3_star(1);
    y                = 1000*(initpoint(2)-e3_star(2));
    X1(ind_T)        = x;
    Y1(ind_T)        = y;
    [thcrit1(ind_T),rcrit1(ind_T)]  = cart2pol(x,y);
    disp([ind_T,length(Tinit1)])
end

Rcrit2  = NaN(size(Tinit2));
thcrit2 = NaN(size(Tinit2));
rcrit2  = NaN(size(Tinit2));

for ind2_T = 1:length(Tinit2)
    initpoint        = deval(BVPper_start,Tinit2(ind2_T));
    GG               = @(RR)tipthre(RR,del,initpoint,Tintg);
    Rinit            = fzero(GG,Rinit);
    Rcrit2(ind2_T)    = Rinit;
    Rinit            = 1.6;
    x                = initpoint(1)-e3_star(1);
    y                = 1000*(initpoint(2)-e3_star(2));
    X2(ind2_T)        = x;
    Y2(ind2_T)        = y;
    [thcrit2(ind2_T),rcrit2(ind2_T)]  = cart2pol(x,y);
    disp([ind2_T,length(Tinit2)])
end

grayscal = 1:-0.01:0;
map = [grayscal',grayscal',grayscal'];
figure;hold on
PCOLOR = pcolor(Rscan-R_star,thper,tip);
PCOLOR.LineStyle = 'none';
PCOLOR.FaceAlpha = 0.3;
colormap(map)
caxis([0 length(ss)]);
plot(Rcrit1-R_star,thcrit1,Rcrit2-R_star,thcrit2,...
    'Color','k','LineWidth',1.5)
xlim([R_end-R_star R_star-R_star])
ylim([-pi pi])
axis([R_end-R_star R_star-R_star -pi pi])
xlabel('$R_- - R_+$','Rotation',0)
ylabel('$\varphi$','Rotation',0)
set(gca,'FontSize',14)
box on
%% phase diagram 
 load('phas_tip_RM_data1.mat')
mult = 1e3;

th_fun1 =  @(RR)interp1(Rcrit1,thcrit1,RR);
th_fun2 =  @(RR)interp1(Rcrit2,thcrit2,RR);

r_fun1  =  @(RR)interp1(Rcrit1,rcrit1,RR);
r_fun2  =  @(RR)interp1(Rcrit2,rcrit2,RR);

R_Test   = 1.8; % threshold for R = 2.5 [2.07];
del = 2.2;

for  n = 1:length(R_Test)
    R_test = R_Test(n);
    par = [R_test,del];
    odefun  = @(t,var)RModefun(var,par);
    
    e0    = [0,0];
    e1   = [R_test/C, 0];
    e2    = [mu,0];
    e3(1) = (del*beta)/((gamma*alpha) - del);
    e3(2) = (R_test-(C*e3(1)))*(beta+e3(1))*(e3(1)-mu)...
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
    
    % intersection points 
    r1 = r_fun1(R_test);
    r2 = r_fun2(R_test);
    
    th1 = th_fun1(R_test);
    th2 = th_fun2(R_test);
    
    N_int1(n) = r1*cos(th1);
    P_int1(n) = r1*sin(th1);
    
    N_int2(n) = r2*cos(th2);
    P_int2(n) = r2*sin(th2);
    
    figure;
    hold on
    box on
%     plot(linspace(3.3,14),10.95*ones(100),'-k')
    plot(...
        e0(1),mult*e0(2),'.',...
        'MarkerSize',20,'Color',[.8 .8 .8])
    plot(...
        e1(1),mult*e1(2),'+',...
        'MarkerSize',7,'LineWidth',1.5,'Color',[.8 .8 .8])
    plot(...
        e3_star(1),mult*e3_star(2),'.',...
        'MarkerSize',20,'Color',[.47 .67 .19])
    
%     plot(...
%         varman1(:,1),mult*varman1(:,2),'k--',...
%         varman2(:,1),mult*varman2(:,2),'k--',...
       plot( varman3(:,1),mult*varman3(:,2),'k--',...
        'LineWidth',1)
    
    plot(varper(1,:),mult*varper(2,:),...
        'LineWidth',1.5,'Color',[.47 .67 .19])
    
    plot(...
        N_int1(n)+e3_star(1),(P_int1(n)+(mult*e3_star(2))),'.',...
        N_int2(n)+e3_star(1),(P_int2(n)+(mult*e3_star(2))),'.',...
        'MarkerSize',15,'Color',[.30 .75 .93]);
    
    axis([-0.5 14 -.001*mult mult*0.025])
    set(gca,'FontSize',14)
    xlabel('$R$')
    ylabel('$P\times10^3$')
end

%% functions

% the flip function 
function [crit] = tipthre(R,del,initpoint,T)
par = [R,del];
odefun=@(t,var)RModefun(var,par);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~,varshot] = ode45(odefun,[0,T],initpoint,opts);
if norm(varshot(end,:)) < 1e-3
    crit = 1;
else
    crit  = -1;
end
end

%the ode fuction of RM model
function [dvar] = RModefun(var,par)
%parameters
% global C gamma beta alpha mu nu delta

RR    = par(1);
delta = par(2);

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
