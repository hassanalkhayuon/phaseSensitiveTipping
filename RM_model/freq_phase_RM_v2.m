% frequency of tipping in terms of phases \varphi
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
% rng(4)
%% time sieries
Tend      = 100;
RR        = 1;  %avreage length of Type-L/H period
PnBin     = .2;
mult      = 1000;
initcond  = [11,0.004];

K = 1;

t_tip   = NaN(1,K);
R_aft   = NaN(1,K);
R_bef   = NaN(1,K);
var_bef = NaN(2,K);

sumbef =0;
numbef =0;

sumaft =0;
numaft =0;

crit = 0;

ind_sim = 0;
% nnn=55;
while (ind_sim < K)
    nnn=nnn+1;
    rng(nnn);
    vars_cl = {'R_c_var','Rcli','R_c','swich','time','cli1','cli','t',...
        'var','NORM'};
    clear(vars_cl{:})
    
    [time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
    climate  = @(tt)interp1(time,cli1,tt);
    R_c      = Rmean + cli .* (Rstar-Rend)/2;
    R        = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;
    n = 1;
    for ind =2: length(R_c)
        Rcli(n) =  R_c(ind-1);
        if (R_c(ind)==R_c(ind-1)) == 0
            n=n+1;
        end
    end
    
    odefun   = @(t,var)RModefun_cli(t,var,R);
    [t,var]  = ode45(odefun,[0:0.1:Tend],initcond,opts);
    
    if t(end) < Tend
        ind_sim = ind_sim+1;
        for ind  = 1:length(var)
            R_c_var(ind) = R(t(ind));
        end
        
        t_tip(ind_sim)    = t(end);
        R_aft(ind_sim)    = R(floor(t_tip(ind_sim)));
        ind_R_aft         = find(Rcli==R_aft(ind_sim),1,'last'); 
        if ind_R_aft >1
            ind_R_bef         = ind_R_aft-1;
        else
            ind_R_bef         = ind_R_aft;
        end
        R_bef(ind_sim)    = Rcli(ind_R_bef);
        
        ind_t_tip_bef    = find(...
            abs(R_c_var - R_bef(ind_sim))< 1e-3,1,'last'); %(Rstar-Rend)/5
        if isfinite(ind_t_tip_bef)
            var_bef(:,ind_sim)     = var(ind_t_tip_bef,:);
        end
        
        % cleaning data 1
        if R_bef(ind_sim) < R_aft(ind_sim)
            odefun1 = @(t,var)RModefun_cli(t,var,@(t)R_bef(ind_sim));
            [t1,var1]  = ode45(odefun1,[0 100],var_bef(:,ind_sim),opts);
            if norm(var1(end,:)) < 1e-3
                R_aft(ind_sim)         = R_bef(ind_sim);
                ind_R_aft              = ind_R_aft-1;
                ind_R_bef              = ind_R_aft-1;
                R_bef(ind_sim)         = Rcli(ind_R_bef);
                ind_t_tip_bef          = find(...
            abs(R_c_var - R_bef(ind_sim))< (Rstar-Rend)/1e-3 ,1,'last'); %(Rstar-Rend)/5 
                if isfinite(ind_t_tip_bef)
                    var_bef(:,ind_sim)     = var(ind_t_tip_bef,:);
                end
            end
            crit = crit+1;
        end
        
        if isfinite(R_aft(ind_sim))
            sumbef = sumbef + R_bef(ind_sim);
            numbef = numbef + 1;
            sumaft = sumaft + R_aft(ind_sim);
            numaft = numaft + 1;
            eref(1) = (del*beta)/((gamma*alpha) - del);
            eref(2) = (R_aft(ind_sim)-(C*eref(1)))*(beta+eref(1))*(eref(1)-mu)/(alpha*(nu + eref(1)));
        else
            eref = NaN(1,2);
        end

        
        [phi_bef(ind_sim),r_bef(ind_sim)]  = cart2pol(...
            var_bef(1,ind_sim)-eref(1),...
            (mult.*(var_bef(2,ind_sim)-eref(2))));
        disp([ind_sim,K])
    end
    
    % cleaning data 2
    if ind_sim ~= 0
        odefun2 = @(t,var)RModefun_cli(t,var,@(t)R_aft(ind_sim));
        [t2,var2]  = ode45(odefun2,[0 100],var_bef(:,ind_sim),opts);
        if (norm (var2(end,:))>1e-3)
            ind_sim = ind_sim -1;
        end
    end
    
    % cleaning data 2.1
    if ind_sim ~= 0
        odefun2 = @(t,var)RModefun_cli(t,var,@(t)R_bef(ind_sim));
        [t2,var2]  = ode45(odefun2,[0 100],var_bef(:,ind_sim),opts);
        if (norm (var2(end,:))<=1e-1)
            ind_sim = ind_sim -1;
        end
    end
     
    % cleaning data 3
    if ind_sim ~= 0
        if norm (var_bef(:,ind_sim))<= 1e-1
            ind_sim = ind_sim -1;
        end
    end  
end
avrR_bef =  sumbef/numbef;
avrR_aft =  sumaft/numaft;
disp([t_tip, max(swich), nnn]);
%% plotting
% load('freq_phase_RM_v4_homo.mat')
% n=1;
% for ind_sim = 1:K
%     if isfinite(phi_bef(ind_sim))
%        t_tip1(n)=t_tip(ind_sim);
%        phi_bef1(n)= phi_bef(ind_sim);
%        n=n+1;
%     end
% end
figure;
subplot(2,2,1)
hold on

plot(R_bef,'r.','MarkerSize',7)
plot([0 K],avrR_bef*[1 1],'-r','LineWidth',2)
plot(R_aft,'b.','MarkerSize',7)
plot([0 K],avrR_aft*[1 1],'-b','LineWidth',2)
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
plot(var_bef(1,:),mult*var_bef(2,:),'k.','Markersize',7)
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
histogram(t_tip,50,'Normalization','probability','FaceColor',[.8 .8 .8],'LineWidth',1.5)

xticks([0 1000 2000 3000 4000 5000])
xticklabels([0 1 2 ])
yticks([0 .05  0.1])
yticklabels([0  .05  0.1])
axis([0 Tend 0 .1])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:100:Tend;
hA.YAxis.MinorTickValues = 0:0.005:.1;
xlabel('tipping time $\times 10^3$')
set(gca,'FontSize',15)
box on

subplot(2,2,4)
hold on
histogram(phi_bef,10,'Normalization','probability','FaceColor',[.8 .8 .8],'LineWidth',1.5)

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi'  '-\pi/2' '0' '\pi/2' ' '})
yticks([0 .1 .2 .3])
yticklabels([0 .1 .2 .3])
axis([0 pi 0 .3])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:0.025*pi:pi;
hA.YAxis.MinorTickValues = 0:0.01:0.5;
xlabel('$\varphi$')
set(gca,'FontSize',15)
box on
%%
% climate variation
function [time,swich,climate,climate1] = Climate(Tend,RR,PnBin)
%climatswich generats a sequence of negative binomial random variables that
%add up to Tend, has propablility of success in singil tril PnBin and
% corresponding number of successes RR
swichtemp = nbininv(rand(1,2*Tend),RR,PnBin);
indend = 1;
while (sum(swichtemp(1:indend))<Tend)
    indend = indend+1;
end
swichtemp1 = swichtemp(1:indend);
ind_sw = 1;
for ind_sw1=1:length(swichtemp1)
    if swichtemp1(ind_sw1)~=0
        swich(ind_sw) = swichtemp1(ind_sw1);
        ind_sw = ind_sw+1;
    end
end
climate1 = NaN(1,sum(swich));
ind_cl = 1;
for ind_sw = 1:length(swich)
    climatampl = rand;
    if rand<=0.5 % good years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = climatampl;
            ind_cl = ind_cl + 1;
        end
    else % bad years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = -climatampl;
            ind_cl = ind_cl + 1;
        end
    end
end
time = 0:0.01:sum(swich);
climate = NaN(size(time));
for ind_cl=1:length(time)-1
    climate(ind_cl) = climate1(floor(time(ind_cl)+1));
end
end

%the ode fuction of RM model
function [dvar] = RModefun_cli(t,var,RRR)
%parameters

R          = RRR(t);
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