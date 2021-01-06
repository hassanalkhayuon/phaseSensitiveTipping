% to simulate the frequency of tipping in terms of phases \varphi when the
% shift pass the homoclinic orbit in RM model

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar  = 2.7;
Rend   = 1.6;
Rmean  = (Rstar+Rend)/2;
del    = 2.2;       %/yr.            Predator death rate in absent of prey
C      = 0.19;      %/yr.            Nonlinear death rate of prey
gamma  = 0.004;     %pred/prey.      Prey-predator conversion rate
beta   = 1.5;       % prey/ha        Predator half-saturating constant
alpha  = 800;       % prey/(pred.yr) Predator saturating kill rate
mu     = 0.03;      % NA             Allee parameter
nu     = 0.003;     % NA             Allee parameter

% opts = odeset('RelTol',1e-3,'AbsTol',1e-4,'Events', @myEvent);
opts = odeset('RelTol',1e-5,'AbsTol',1e-10,'Events', @myEvent);
rng(4)
%% Simulations

Tend      = 5000;
RR        = 1;  %avreage length of Type-L/H period
PnBin     = .2;
mult      = 1000;


K = 1000;

t_tip   = NaN(1,K);
R_aft   = NaN(1,K);
R_bef   = NaN(1,K);
var_bef = NaN(2,K);

ind_sim  = 1;
ind_rate = 1;
ind_bif  = 1;
count = 0;
while (ind_sim < K)
    
    vars_cl = {'T','ind_T','t','var','R'};
    clear(vars_cl{:})
    
    T(1)           = nbininv(rand,RR,PnBin);
    R(1)           = NaN;
    initcond(:,1)  = [11;0.004];
    while (T(1) == 0)
        T(1) = nbininv(rand,RR,PnBin);
    end
    
    ind_T = 2;
    
    NORM = inf;
    Ttip = 0;
    while and(sum(T)<Tend , NORM>1e-3)
        T(ind_T) = nbininv(rand,RR,PnBin);
        while (T(ind_T) == 0)
            T(ind_T) = nbininv(rand,RR,PnBin);
        end
        R(ind_T) = Rend + rand*(Rstar - Rend);
        odefun   = @(t,var)RModefun(var,R(ind_T));
        tspan    = [sum(T(1:ind_T - 1)),sum(T(1:ind_T))];
        [t,var]  = ode45(odefun,tspan,initcond(:,ind_T-1),opts);
        initcond(:,ind_T) = var(end,:);
        NORM = norm(initcond(:,ind_T));
        Ttip = Ttip + T(ind_T);
        ind_T = ind_T + 1;
    end
    
    % simulations outcomes
    R_bef(ind_sim)      = R(ind_T - 2);
    R_aft(ind_sim)      = R(ind_T - 1);
    t_tip(ind_sim)      = Ttip;
    
    eref(1) = (del*beta)/((gamma*alpha) - del);
    eref(2) = (R_aft(ind_sim)-(C*eref(1)))*(beta+eref(1))*(eref(1)-mu)/(alpha*(nu + eref(1)));
    
    var_bef(:,ind_sim)  = initcond(:,ind_T - 2);
    
    [phi_bef(ind_sim),r_bef(ind_sim)]  = cart2pol(...
        var_bef(1,ind_sim)-eref(1),...
        (mult.*(var_bef(2,ind_sim)-eref(2))));
    
    % Testing data 1
    odefun_bef   = @(t,var)RModefun(var,R_bef(ind_sim));
    odefun_aft   = @(t,var)RModefun(var,R_aft(ind_sim));
    [~,var_b]  = ode45(odefun_bef,[0,10],var_bef(:,ind_sim),opts);
    [~,var_a]  = ode45(odefun_aft,[0,10],var_bef(:,ind_sim),opts);
    NORM_bef(ind_sim) = norm(var_b(end,:)); %NORM_bef >=0.1;
    NORM_aft(ind_sim) = norm(var_a(end,:)); %NORM_aft <=1e-3;
    if and(NORM_bef(ind_sim) > 1e-1,NORM_aft(ind_sim)<1e-3)
        if (R_aft(ind_sim)>=2.57)            
            var_bef_bif(:,ind_bif) = var_bef(:,ind_sim);
            
            [phi_bef_bif(ind_bif),r_bef_bef(ind_bif)]  = cart2pol(...
                var_bef_bif(1,ind_bif)-eref(1),...
                (mult.*(var_bef_bif(2,ind_bif)-eref(2))));
            
            t_tip_bif (ind_bif) = Ttip;
            
            ind_bif = ind_bif + 1;   
        else
            if  and((R_aft(ind_sim) < R_bef(ind_sim)),...
                    and(norm (var_bef(:,ind_sim))> 0.07, R_aft(ind_sim) < 2.65) )
                
                var_bef_rate(:,ind_rate) = var_bef(:,ind_sim);
                
                [phi_bef_rate(ind_rate),r_bef_rate(ind_rate)]  = cart2pol(...
                    var_bef_rate(1,ind_rate)-eref(1),...
                    (mult.*(var_bef_rate(2,ind_rate)-eref(2))));
                
                t_tip_rate (ind_rate) = Ttip;
                
                ind_rate = ind_rate + 1;
            else
                ind_sim = ind_sim -1;
                count = count +1;
            end
            
        end
        ind_sim=ind_sim+1;
    end
    
    disp([ind_sim,K])
end

%% plotting

figure;
subplot(2,2,1)
hold on

plot(R_bef,'r.','MarkerSize',7)
% plot([0 K],mean(R_bef(1:end))*[1 1],'-r','LineWidth',2)
plot(R_aft,'b.','MarkerSize',7)
% plot([0 K],mean(R_aft(1:end))*[1 1],'-b','LineWidth',2)
xticks([0 200 400 600 800 1000])
xticklabels([0 200 400 600])
yticks([0 .5 1 1.5 2 2.5 3 3.5 4])
yticklabels([0 .5 1 1.5 2 2.5 3 3.5 4])
axis([0 K 1.5 3])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:20:1000;
hA.YAxis.MinorTickValues = 0:0.05:5;
set(gca,'FontSize',15)
box on
xlabel('Simulation index')

subplot(2,2,2)
hold on
plot(var_bef_rate(1,:),mult*var_bef_rate(2,:),'k.','Markersize',7)
plot(var_bef_bif(1,:),mult*var_bef_bif(2,:),'.','Markersize',7,'Color',[.7 .7 .7])


axis([0 15 0 35])

xticks([0 2 4 6 8 10])
xticklabels([0 2 4 6 8])
yticks([0 5 10 15 20 25 30 35])
yticklabels([0 5 10 15 20 25 30])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:0.5:10;
hA.YAxis.MinorTickValues = 0:1:35;
xlabel('$N$')
ylabel('\tilde{P}','Rotation',0)
set(gca,'FontSize',15)
box on

subplot(2,2,3)
hold on
histogram(...
    t_tip_bif,60,'Normalization','probability',...
    'FaceColor',[.8 .8 .8],'LineWidth',1.5)

histogram(...
    t_tip,100,'Normalization','probability',...
    'FaceColor','k','LineWidth',1.5)

xticks([0 1000 2000 3000 4000 5000])
xticklabels([0 1 2 ])
yticks([0 .05  0.1])
yticklabels([0  .05  0.1])
axis([0 1000 0 .1])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:100:Tend;
hA.YAxis.MinorTickValues = 0:0.005:.1;
xlabel('tipping time $\times 10^3$')
set(gca,'FontSize',15)
box on

subplot(2,2,4)
hold on

% histogram(...
%     phi_bef,60,'Normalization','probability',...
%     'FaceColor',[.8 .8 .8],'LineWidth',1.5)


histogram(...
    phi_bef_bif,60,'Normalization','count',...
    'FaceColor',[.8 .8 .8],'LineWidth',1.5)
histogram(...
    phi_bef_rate,60,'Normalization','count',...
    'FaceColor','k','LineWidth',1.5)

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi'  '-\pi/2' '0' '\pi/2' ' '})
yticks([0 50 100 150 200])
yticklabels({0 ''  0.1 '' 0.2})
axis([-pi pi 0 250])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:0.025*pi:pi;
hA.YAxis.MinorTickValues = 0:0.01:0.5;
xlabel('$\varphi$')
set(gca,'FontSize',15)
box on
%%

%the ode fuction of RM model
function [dvar] = RModefun(var,R)
%parameters

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

% myEvent to stop the integration when it tips
function [value, isterminal, direction] = myEvent(t,y)
TOL = 1e-3;
value      = norm(y)<TOL;
isterminal = 1;   % Stop the integration
direction  = 0;   % approch zero from either diractions?
end