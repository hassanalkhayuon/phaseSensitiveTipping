
warning off
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');



% parameters:

Rlow    = 2.3;       % /yr.              Prey intrinsic growth rate
Rupp    = 3.3;
Rmid    = 3;
Rend    = Rupp;
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

T        = 100;
 
 opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
%% step 1 defining periodic orbits

% periodic orbit for R = Rlow
mult     = 1000;
initcond = [10   0.001];
Tper     =  9;

ivpfun   = @(t,var)Mayodefun(var,Rlow);
perfun   = @(t,var,Tper)Tper*Mayodefun(var,Rlow);

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

per_low = bvp5c(perfun,BC,solinit);

% periodic orbit for R=Rmid
mult     = 1000;
initcond = [10   0.001];
Tper     =  8.6;

ivpfun   = @(t,var)Mayodefun(var,Rmid);
perfun   = @(t,var,Tper)Tper*Mayodefun(var,Rmid);

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

per_mid = bvp5c(perfun,BC,solinit);

% periodic orbit for R=Rupp
load('data\perinitstart_may.mat')

per_upp = BVPper_start;

N_initconds = 1000; % number of initial conditions on the limit cycle
ss = linspace(0,1,N_initconds);

varlow = deval(per_low,ss);
varmid = deval(per_mid,ss);
varupp = deval(per_upp,ss);

figure; hold on
plot(...
    varlow(1,:)+1,mult*varlow(2,:)+6,...
    'LineWidth',3,'Color',[0.47 0.67 0.19],'Marker','.')
% plot(...
%     varmid(1,:),mult*varmid(2,:),...
%     'LineWidth',2,'Color','g','Marker','.')
plot(...
    varupp(1,:),mult*varupp(2,:),...
    'LineWidth',2,'Color','r','Marker','.')
%%
Rlow       = 2.3; 
ivpfun   = @(t,var)Mayodefun(var,Rlow);
%[~,var]  = ode45(ivpfun,[0 5],[10.5 0.0001],opts);
[~,var]  = ode45(ivpfun,[0 7], [1.1    0.004] ,opts);
%figure
plot(...
    var(7:end,1)+1,mult*var(7:end,2)+6,'--k',...
    'LineWidth',1.5)
%%
point1 = [0.095 4.946];
point2 = [-1.4 4];
point3 = [-3 0];
yy = @(xx)spline([point1(1), point2(1), point3(1)],[point1(2), point2(2), point3(2)],xx);
figure(2);
fplot(yy,[-3 0.095])
%% replicating the two periodic orbits and \theta(R)

% the start periodic orbit
Nupp_top = [(varupp(1,684:1000)),(varupp(1,2:276))];
Pupp_top = mult*[(varupp(2,684:1000)),(varupp(2,2:276))];

Nupp_dow = varupp(1,276:684);
Pupp_dow = mult*varupp(2,276:684);

Pstar_top = @(NN)interp1(Nupp_top,Pupp_top,NN);
Pstar_dow = @(NN)interp1(Nupp_dow,Pupp_dow,NN);

% the middle periodic orbit 
% 681  228
Nmid_top = [(varmid(1,681:1000)),(varmid(1,2:228))];
Pmid_top = mult*[(varmid(2,681:1000)),(varmid(2,2:228))];

Nmid_dow = varmid(1,228:681);
Pmid_dow = mult*varmid(2,228:681);

Pmidl_top = @(NN)interp1(Nmid_top,Pmid_top,NN);
Pmidl_dow = @(NN)interp1(Nmid_dow,Pmid_dow,NN);

% the end periodic orbit
Nlow_top = [varlow(1,871:999),varlow(1,1:214)];
Plow_top = mult*[varlow(2,871:999),varlow(2,1:214)];

Nlow_dow = varlow(1,214:871);
Plow_dow = mult*varlow(2,214:871);

Pend_top = @(NN)interp1(Nlow_top,Plow_top,NN);
Pend_dow = @(NN)interp1(Nlow_dow,Plow_dow,NN);

% threshold 
odefun  = @(t,var)Mayodefun(var,Rend);
[N_eq] = May_eq(Rend,q);
e5(1) =  N_eq(3);
e5(2) = (N_eq(3) + eps)/q;
G   = @(var)Mayodefun(var,Rend);
JJ = MyJacobian(G,e5);
[eigvic,eigval]=eig(JJ);
if eigval(2,2)<0
    pert = eigvic(:,2);
else
    pert = eigvic(:,1);
end
pert = 1e-4*pert';
maninit = e5 + pert;
[~,varman]   = ode45(odefun,[10,0],maninit,opts);
Pthe = @(NN)interp1(varman(:,1),mult*varman(:,2),NN);

% % figure; hold on
% plot(...
%     varlow(1,:),mult*varlow(2,:),...
%     'LineWidth',2,'Color','b','Marker','.')
% plot(...
%     Nlow_top,Plow_top,'o')
%% scan in (N,P)-plane

grid_N = 1000;
grid_P = 2*grid_N;
N_scan = linspace(0,15,grid_N);
P_scan = linspace(0,35,grid_P);
colormat = zeros(grid_P,grid_N);

for ind_N = 1:grid_N
    N  = N_scan(ind_N);

    
    % Region 2.1, green color
    if and(N>=0.4218,N<=13.26)
        Pmax2_1 = min(Pthe(N),Pstar_top(N));
        Pmin2_1 = max(Pend_top(N),min(Pstar_dow(N),Pmidl_dow(N)));
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_1, P>Pmin2_1)
                colormat(ind_P,ind_N) = 2;
            end
        end
    end
    
    % Region 2.2
    if and(N>=3.553,N<6.066)
        Pmax2_2 = Pend_dow(N);
        Pmin2_2 = min(Pstar_dow(N),Pmidl_dow(N));
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_2, P>Pmin2_2)
                colormat(ind_P,ind_N) = 2;
            end
        end
    end
%     
    % Region 2.3
    if N<=3.553
        Pmax2_3 = Pstar_dow(N);
        Pmin2_3 = min(Pend_dow(N),Pmidl_dow(N));
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_3, P>Pmin2_3)
                colormat(ind_P,ind_N) = 2;
            end
        end
    end
%     
%     % Region 2.4
    if N<0.21
        Pmax2_4 = max(Pend_top(N),Pmidl_top(N));
        Pmin2_4 = Pend_dow(N);
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_4, P>Pmin2_4)
                colormat(ind_P,ind_N) = 2;
            end
        end
    end
%     
    % Region 2.5
    if N>=0.21
        Pmax2_5 = min(max(Pend_top(N),Pmidl_top(N)),Pthe(N));
        Pmin2_5 = Pstar_top(N);
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax2_5, P>Pmin2_5)
                colormat(ind_P,ind_N) = 2;
            end
        end
    end
    
        % Region 1, red color 
    if N<13.26
        Pmax1_1 = max(Pstar_top(N),Pmidl_top(N));
        Pmin1_1 = Pthe(N);%max(max(Pend_top(N),Pthe(N)),Pstar_dow(N));
        parfor ind_P = 1:grid_P
            P = P_scan(ind_P);
            if and(P<=Pmax1_1, P>Pmin1_1)
                colormat(ind_P,ind_N) = 1;
            end
        end
    end
    
%     % Region 1.2
%     if N<13.26
%         Pmax1_2 = Pend_top(N);
%         Pmin1_2 = max(Pthe(N),Pstar_top(N));
%         parfor ind_P = 1:grid_P
%             P = P_scan(ind_P);
%             if and(P<=Pmax1_2, P>Pmin1_2)
%                 colormat(ind_P,ind_N) = 1;
%             end
%         end
%     end
end

%% plotting 

figure;
map = [1,1,1;1 0 0;.47,.67,.19];
PCOLOR = pcolor(N_scan,P_scan,colormat);
PCOLOR.FaceAlpha = 0.1;
PCOLOR.LineStyle = 'none';
colormap(map)
hold on

fplot(Pstar_top,[0 15])
fplot(Pstar_dow,[0 15])

fplot(Pmidl_top,[0 15])
fplot(Pmidl_dow,[0 15])

fplot(Pend_top,[0 15])
fplot(Pend_dow,[0 15])

fplot(Pthe,[0,15])

caxis([0 2]);
%% functions 

function [dvar] = Mayodefun(var,R)

%parameters

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

