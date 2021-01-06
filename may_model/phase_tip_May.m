% To find the tipping phases of the may model.
% the starting point is at R = 3.3, q = 205. q is going to be fixes through
% out and R is going to be decreased up to R_H = 1.67.

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');
% parameters:

R_star  = 2.8;       % /yr.              Prey intrinsic growth rate
R_end   = 1.65;
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%% initiation

N_initconds = 500; % number of initial conditions on the limit cycle
ss = linspace(0,1,N_initconds);
load('data\perinitstart_may.mat')
varper = deval(BVPper_start,ss);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

[N_eq] = May_eq(R_star,q);
e3_star(1) =  N_eq(1);
e3_star(2) = (N_eq(1) + eps)/q;

[thper,rper]  = cart2pol(varper(1,:)-e3_star(1),(1000*(varper(2,:)-e3_star(2))));
AA = sortrows([thper',rper']);
thper = AA(:,1);
rper  = AA(:,2);
Rscan = linspace(R_end,R_star,1000);

Tint  = 100;
e0    = [0,eps/q];
figure; hold on

plot(...
    rper(368:483).*cos(thper(368:483))+e3_star(1),...
    rper(368:483).*sin(thper(368:483))+1000*e3_star(2),...
    'r-','LineWidth',2)
plot(...
    rper(1:483).*cos(thper(1:483))+e3_star(1),...
    rper(1:483).*sin(thper(1:483))+1000*e3_star(2),...
    'r-','LineWidth',2)
%  plot([2.069 2.006],[21.68 21.33],'r-','LineWidth',2)
% 
% plot(0,0,'.k','MarkerSize',15)
% plot(e3_star(1),e3_star(2),'k.','MarkerSize',15)
% plot(varper(1,:),varper(2,:),'ko','LineWidth',2)
% plot(varper(1,497)-e3_star(1),varper(2,497)-e3_star(2),'r.','MarkerSize',15)
% plot([varper(1,497)-e3_star(1) 0],[varper(2,497)-e3_star(2) 0],'--k')
% plot([0 10],[0 0],'k--')

%% scan 
tic
for indR = 1:length(Rscan)
    R = Rscan(indR);
    TIPind = 1;
    parfor indth = 1:length(thper)
        N = rper(indth)*cos(thper(indth))+ e3_star(1);
        P = (rper(indth)*sin(thper(indth))/1000)+ e3_star(2);
        init    = [N;P]
        par     = [R,q];
        odefun  = @(t,var)Mayodefun(var,par);
        [t,var] = ode45(odefun,[0 Tint],init,opts);
        if norm(var(end,:) - e0)<=1e-4
            tip(indth,indR) = 1;
            TIPind = TIPind+1;
        else
            tip(indth,indR) = 0;
        end
    end
    tip(:,indR)= tip(:,indR) * TIPind;
    disp([indR,length(Rscan)])
end

% grayscal = 1:-0.01:0;
% map = [grayscal',grayscal',grayscal'];
% figure;hold on
% PCOLOR = pcolor(Rscan-R_star,thper,tip);
% PCOLOR.LineStyle = 'none';
% PCOLOR.FaceAlpha = 0.3;
% colormap(map)
% caxis([0 10]);


%% threshold
Ttop    = 1;
Tend1   = 0.239523952395240; % find Tend1 for R = 1.67 
Tend2   = 0.458345834583458;
Tinit1  = Tend2:0.0001:1;
Tinit2  = 0:0.0001:Tend1;
Tintg   = 100;
Rinit   = R_end;
Rcrit1  = NaN(size(Tinit1));
thcrit1 = NaN(size(Tinit1));
rcrit1  = NaN(size(Tinit1));
for ind_T = 1:length(Tinit1)
    initpoint        = deval(BVPper_start,Tinit1(ind_T));
    GG               = @(RR)tipthre(RR,q,initpoint,Tintg);
    Rinit            = fzero(GG,Rinit);
    Rcrit1(ind_T)    = Rinit;
    x                = initpoint(1)-e3_star(1);
    y                = 1000*(initpoint(2)-e3_star(2));
    X1(ind_T)        = x;
    Y1(ind_T)        = y;
    [thcrit1(ind_T),rcrit1(ind_T)]  = cart2pol(x,y);
    disp([ind_T,length(Tinit1)])
end

Rinit   = R_end;
Rcrit2  = NaN(size(Tinit2));
thcrit2 = NaN(size(Tinit2));
rcrit2  = NaN(size(Tinit2));

for ind_T = 1:length(Tinit2)
    initpoint       = deval(BVPper_start,Tinit2(ind_T));
    GG              = @(RR)tipthre(RR,q,initpoint,Tintg);
    Rinit           = fzero(GG,Rinit);
    Rcrit2(ind_T)   = Rinit;
    x               = initpoint(1)-e3_star(1);
    y               = 1000*(initpoint(2)-e3_star(2));
    X2(ind_T)       = x;
    Y2(ind_T)       = y;
    [thcrit2(ind_T),rcrit2(ind_T)]  = cart2pol(x,y);
    disp([ind_T,length(Tinit2)])

end

kk = find(thcrit2<0,1);
thcrit_fun2 =  @(RR)interp1(Rcrit2(1:kk-1),thcrit2(1:kk-1),RR);
thcrit_fun3 =  @(RR)interp1(Rcrit2(kk:end),thcrit2(kk:end),RR);
thcrit_fun1 =  @(RR)interp1(Rcrit1,thcrit1,RR);
%%
figure;
hold on

grayscal = 1:-0.01:0;
map = [grayscal',grayscal',grayscal'];

PCOLOR = pcolor(Rscan-R_star,thper,tip);
PCOLOR.LineStyle = 'none';
PCOLOR.FaceAlpha = 0.3;
colormap(map)
caxis([0 N_initconds]);
Rscan = linspace(R_end,R_star,10000);
plot(...
    Rscan-R_star,thcrit_fun1(Rscan),'k-',...
    Rscan-R_star,thcrit_fun2(Rscan),'k-',...
    Rscan-R_star,thcrit_fun3(Rscan),'k-',...
    'LineWidth',1.5)
% plot(...
%     [3.3-R_star,1.67-R_star],[pi pi],'--k',...
%     [3.3-R_star,1.67-R_star],[-pi -pi],'--k',...
%     'LineWidth',0.5)

axis([-1.63 0 -pi pi])
xlabel('$R_\mathrm{start} - R$','Rotation',0)
ylabel('$\beta$','Rotation',0)
set(gca,'FontSize',14)
box on


%% visual test
load('data\visual_test_intersections.mat');
th_fun1 =  @(RR)interp1(Rcrit1,thcrit1,RR);
th_fun2 =  @(RR)interp1(Rcrit2(1:kk-1),thcrit2(1:kk-1),RR);
th_fun3 =  @(RR)interp1(Rcrit2(kk:end),thcrit2(kk:end),RR);

r_fun1  = @(RR)interp1(Rcrit1,rcrit1,RR);
r_fun2  = @(RR)interp1(Rcrit2(1:kk-1),rcrit2(1:kk-1),RR);
r_fun3 =  @(RR)interp1(Rcrit2(kk:end),rcrit2(kk:end),RR);


R_Test = 3.3;%2.4082:-0.08:R_end;
load('data\visual_test_intersections.mat')
for  n = 1:length(R_Test)
    R_test = R_Test(n);
    par   = [R_test,q];
    odefun  = @(t,var)Mayodefun(var,par);
    
%     equilibria
    
    [N_eq] = May_eq(par(1),par(2));
    
    e0    = [0,eps/par(2)];
    e1    = [par(1)/C, 0];
    e2    = [mu,0];
    
    e3(1) =  N_eq(1);
    e3(2) = (N_eq(1) + eps)/par(2);
    
    e4(1) =  N_eq(2);
    e4(2) = (N_eq(2) + eps)/par(2);
    
    e5(1) =  N_eq(3);
    e5(2) = (N_eq(3) + eps)/par(2);
    
    G   = @(var)Mayodefun(var,par);
    JJ = MyJacobian(G,e5);
    
    [eigvic,eigval]=eig(JJ);
    if eigval(2,2)<0
        pert = eigvic(:,2);
    else
        pert = eigvic(:,1);
    end
    pert = 1e-4*pert';
    
    
    JJ1 = MyJacobian(G,e1);
    
    [eigvic,eigval]=eig(JJ1);
    if eigval(2,2)<0
        pert1 = eigvic(:,2);
    else
        pert1 = eigvic(:,1);
    end
    pert1 = 1e-4*pert1';
    
    maninit1 = e5 + pert;
    maninit2 = e5 - pert;
    maninit3 = e1 + pert1;
    maninit4 = e1 - pert1;
    
    [~,varman1]   = ode45(odefun,[20,0],maninit1,opts);
    [~,varman2]   = ode45(odefun,[20,0],maninit2,opts);
    [~,varman3]   = ode45(odefun,[0,50],maninit3,opts);
    [~,varman4]   = ode45(odefun,[0,50],maninit4,opts);
    
    % intersection points 
    
    if R_test <= Rcrit2(kk)
        r1 = r_fun1(R_test);
        r2 = r_fun3(R_test);
        
        th1 = th_fun1(R_test);
        th2 = th_fun3(R_test);
        
        N_int1(n) = r1*cos(th1);
        P_int1(n) = r1*sin(th1);
        
        N_int2(n) = r2*cos(th2);
        P_int2(n) = r2*sin(th2);
    else
        r1 = r_fun1(R_test);
        r2 = r_fun2(R_test);
        
        th1 = th_fun1(R_test);
        th2 = th_fun2(R_test);
        
        N_int1(n) = r1*cos(th1);
        P_int1(n) = r1*sin(th1);
        
        N_int2(n) = r2*cos(th2);
        P_int2(n) = r2*sin(th2);
    end

figure
%     cla
    hold on
    plot(varper(1,:),1000*varper(2,:),'-b','LineWidth',2)
        plot(...
            e3_star(1),1000*e3_star(2),'.',...
            'MarkerSize',20,'Color',[.47 .67 .19])
        plot(...
            e0(1),e0(2),'ok',...
            'MarkerSize',4,'LineWidth',3.5)
        plot(...
            e2(1),e2(2),'+k',...
            e5(1),e5(2),'+k',...
            'MarkerSize',8,'LineWidth',1)
        plot(...
            e1(1),e1(2),'ok',...
            e3(1),e3(2),'ok',...
            e4(1),e4(2),'ok',...
            'MarkerSize',6,'LineWidth',1)
    plot(...
        varman1(:,1),1000*varman1(:,2),'--k',...
        varman2(:,1),1000*varman2(:,2),'--k')
        
    plot(...
        varman3(:,1),1000*varman3(:,2),'--k',...
        varman4(:,1),1000*varman4(:,2),'--k')
    
    plot(...
        N_int1(n)+e3_star(1),1000*(P_int1(n)+e3_star(2)),'.',...
        N_int2(n)+e3_star(1),1000*(P_int2(n)+e3_star(2)),'.',...
        'MarkerSize',15,'Color',[.30 .75 .93]);
    
    axis([-1 15 -2 40])
    set(gca,'FontSize',14)
    xlabel('$N$')
    ylabel('$P\times10^3$')

    box on
%     grid on
%     grid minor
%     axis([0 15 0 15])
    
%     subplot(2,1,2)
%     cla
%     hold on
%     plot(...
%     Rscan,thcrit_fun1(Rscan),'k-',...
%     Rscan,thcrit_fun2(Rscan),'k-',...
%     Rscan,thcrit_fun3(Rscan),'k-',...
%     'LineWidth',1.5)
% 
%     plot(...
%         R_test*[1 1], [-4 4],'--','Color',[.7 .7 .7])


%     drawnow
%     pause(0.25)
end

%% typical trajectories

initcond = [3.227,19.82/1000];
par      = [3.3,q];
odefun   = @(t,var)Mayodefun(var,par);
[~,vartemp] = ode45(odefun,[0,20],initcond ,opts);
plot(vartemp(:,1),1000*vartemp(:,2),'k-','LineWidth',1.5)

%% functions

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

% the ode functiom
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