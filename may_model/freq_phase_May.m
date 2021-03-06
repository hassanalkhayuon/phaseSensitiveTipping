% frequency of tipping in terms of phases \varphi
warning off
% clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar   = 3.3;       % /yr.              Prey intrinsic growth rate
Rend    = 2;
Rmean   = (Rstar+Rend)/2;
q       = 205;       % prey/pred         Minimum prey biomass.
C       = 1.75/8;    % /yr.              Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr               Rate of predator population.
mu      = 0.03;      % NA                Allee parameter.
nu      = 0.003;     % NA                Allee parameter.
eps     = 0.031;     % NA                Artificial parameter.

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
rng(4)   % random seeding 
%% time sieries
Tend      = 400;
RR        = 10;  %avreage length of Type-L/H period
PnBin     = .5;
mult      = 1000;
initcond  = [11,0.004];

K = 10;

t_tip   = NaN(1,K);
R_aft   = NaN(1,K);
R_bef   = NaN(1,K);
var_bef = NaN(2,K);

sumbef =0;
numbef =0;

sumaft =0;
numaft =0;

for ind_sim = 1:K
    vars_cl = {'R_c_var','Rcli','R_c','swich','time','cli1','cli','t',...
        'var','NORM'};
    clear(vars_cl{:})
    
    [time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
    climate  = @(tt)interp1(0:1:sum(swich)-1,cli,tt);
    R_c      = Rmean + cli .* (Rstar-Rend)/2;
    R        = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;
    n = 1;
    for ind =2: length(R_c)
        Rcli(n) =  R_c(ind-1);
        if (R_c(ind)==R_c(ind-1)) == 0
            n=n+1;
        end
    end
    
    odefun   = @(t,var)Mayodefun_cli(t,var,R);
    [t,var]  = ode45(odefun,[0 Tend],initcond,opts);
    Nt = @(tt)interp1(t,var(:,1),tt);
    Pt = @(tt)interp1(t,var(:,2),tt);
    parfor ind  = 1:length(var)
        NORM(ind)    = norm (var(ind,:) - [0 eps/q]);
        R_c_var(ind) = R(t(ind));
    end
    ind_t_tip_aft   = find(NORM<=1e-3,1,'first');
    if isfinite(ind_t_tip_aft)
        t_tip(ind_sim)         = t(ind_t_tip_aft);
        R_aft(ind_sim)         = R(ceil(t_tip(ind_sim)));
        ind_R_aft              = find(Rcli==R_aft(ind_sim));
        ind_R_bef              = ind_R_aft-1;
        R_bef(ind_sim)         = Rcli(ind_R_bef);
        ind_t_tip_bef          = find(R_c_var == R_bef(ind_sim),1,'last');
        if isfinite(ind_t_tip_bef)
            var_bef(:,ind_sim)     = var(ind_t_tip_bef,:);
        end
    end
       
    % cleaning the data 1.
    if R_bef(ind_sim) < R_aft(ind_sim)  
        odefun1 = @(t,var)Mayodefun_cli(t,var,@(t)R_bef(ind_sim));
        [t1,var1]  = ode45(odefun1,[0 100],var_bef(:,ind_sim),opts);
        if norm (var1(end,:) - [0 eps/q])<1e-3
            R_aft(ind_sim)         = R_bef(ind_sim);
            ind_R_aft              = ind_R_aft-1;
            ind_R_bef              = ind_R_aft-1;
            R_bef(ind_sim)         = Rcli(ind_R_bef);
            ind_t_tip_bef          = find(R_c_var == R_bef((ind_sim)),1,'last');
            if isfinite(ind_t_tip_bef)
                var_bef(:,ind_sim)     = var(ind_t_tip_bef,:);
            end
        end
    end
    
    % cleaning the data 2.
    odefun2 = @(t,var)Mayodefun_cli(t,var,@(t)R_aft(ind_sim));
    [t2,var2]  = ode45(odefun2,[0 100],var_bef(:,ind_sim),opts);
    if norm (var2(end,:) - [0 eps/q])>1e-3
        R_aft(ind_sim)         = NaN;
        ind_R_aft              = NaN;
        ind_R_bef              = NaN;
        R_bef(ind_sim)         = NaN;
        ind_t_tip_bef          = NaN;
        var_bef(:,ind_sim)     = NaN;
    end

    if isfinite(R_aft(ind_sim))
        sumbef = sumbef + R_bef(ind_sim);
        numbef = numbef + 1;
        sumaft = sumaft + R_aft(ind_sim);
        numaft = numaft + 1;                      % to calculate the avrage 
        [N_eq] = May_eq(R_aft(ind_sim),q);
        if isfinite(N_eq(1))
            eref(1) = N_eq(1);                    % the equilibria for phi
            eref(2) = (N_eq(1) + eps)/q;
        elseif isfinite(N_eq(2))
            eref(1) = N_eq(2);
            eref(2) = (N_eq(2) + eps)/q;
        end
    end
%     
    % the phase phi
    [phi_bef(ind_sim),r_bef(ind_sim)]  = cart2pol(...
        var_bef(1,ind_sim)-eref(1),...
        (mult.*(var_bef(2,ind_sim)-eref(2))));
    disp([ind_sim,K])
end
avrR_bef =  sumbef/numbef;
avrR_aft =  sumaft/numaft;

%% plotting
load('freq_climate1.mat')
n=1;
for ind_sim = 1:K
    if isfinite(phi_bef(ind_sim))
       t_tip1(n)=t_tip(ind_sim);
       phi_bef1(n)= phi_bef(ind_sim);
       n=n+1;
    end
end
figure;
subplot(2,2,1)
hold on

plot(R_bef,'r.','MarkerSize',7)
plot([0 K],avrR_bef*[1 1],'-r','LineWidth',2)
plot(R_aft,'b.','MarkerSize',7)
plot([0 K],avrR_aft*[1 1],'-b','LineWidth',2)
xticks([0 200 300 500 700 900])
xticklabels([0 100 300 500])
yticks([2 2.5 3 3.5 4])
yticklabels([2 2.5 3 3.5 4])
axis([0 K 1.8 4])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:20:1000;
hA.YAxis.MinorTickValues = 1.5:0.05:5;
set(gca,'FontSize',15)
box on
xlabel('Simulation index')

subplot(2,2,2)
hold on
plot(var_bef(1,:),mult*var_bef(2,:),'k.','Markersize',7)
axis([0 15 0 35])

xticks([0 3 6 9 12 15])
xticklabels([0 3 6 9 12])
yticks([0 5 15 25 35])
yticklabels([0 5 15 25 ])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:0.5:15;
hA.YAxis.MinorTickValues = 0:1:35;
xlabel('$N$')
ylabel('\tilde{P}','Rotation',0)
set(gca,'FontSize',15)
box on

subplot(2,2,3)
hold on
histogram(t_tip1,50,'Normalization','probability','FaceColor',[.8 .8 .8],'LineWidth',1.5)

xticks([0 100 200 300 400])
xticklabels([0 100 200])
yticks([0 .02 .04 .06 .08 .1])
yticklabels([0 .02 .04 .06 .08 .1])
axis([0 400 0 .1])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = 0:10:400;
hA.YAxis.MinorTickValues = 0:0.005:0.1;
xlabel('Time until tipping')
set(gca,'FontSize',15)
box on

subplot(2,2,4)
hold on
histogram(phi_bef1,50,'Normalization','probability','FaceColor',[.8 .8 .8],'LineWidth',1.5)

xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi' '-\pi/2' '0' '\pi/2' ' '})
yticks([0 .1 .2 .3 .4])
yticklabels([0 .1 .2 .3 .4])
axis([-pi pi 0 .4])
hA = gca;
set(gca,'XMinorTick','on','YMinorTick','on');
hA.XAxis.MinorTickValues = -pi:0.05*pi:pi;
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
    if mod(ind_sw,2)==1 % good years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = climatampl;
            ind_cl = ind_cl + 1;
        end
    elseif mod(ind_sw,2)==0 % bad years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = - climatampl;
            ind_cl = ind_cl + 1;
        end
    end
end
time = 0:0.001:sum(swich);
climate = NaN(size(time));
for ind_cl=1:length(time)-1
    climate(ind_cl) = climate1(floor(time(ind_cl)+1));
end
end

% ode of the May Model
function [dvar] = Mayodefun_cli(t,var,RRR)

%parameters
R       = RRR(t);
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