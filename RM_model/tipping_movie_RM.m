% tipping_movie_RM preduces a movie to illustrate phase-sensetiv tipping in
% the RM model

warning off
clc
clear

addpath('C:\Users\halkhayuon\Dropbox\01_researchprojects\03-codes\01_MATLAB')
addpath('C:\Users\ha317\Dropbox\01_researchprojects\03-codes\01_MATLAB')
set(0,'defaulttextInterpreter','latex');

% parameters:

Rstar     = 2.5;
Rend      = 1.6;
del        = 2.2;    %/yr.            Predator death rate in absent of prey
C          = 0.19;   %/yr.            Nonlinear death rate of prey
gamma      = 0.004;  %pred/prey.      Prey-predator conversion rate
beta       = 1.5;    % prey/ha        Predator half-saturating constant
alpha      = 800;    % prey/(pred.yr) Predator saturating kill rate
mu         = 0.03;   % NA             Allee parameter
nu         = 0.003;  % NA             Allee parameter

opts1 = odeset('RelTol',1e-5,'AbsTol',1e-10);
opts = odeset('RelTol',1e-16,'AbsTol',1e-16);
%% 1) climate variability
rng(301)
Tend      = 100;
RR        = 1;  %avreage length of Type-L/H period
PnBin     = .2;
mult      = 1000;

T(1) = nbininv(rand,RR,PnBin);
R(1) = Rend + rand*(Rstar - Rend);
initcond(:,1)  = [3;0.002];
while (T(1) == 0)
    T(1) = nbininv(rand,RR,PnBin);
end

odefun   = @(t,var)RModefun(var,R(1));
tspan    = [0,T(1)];
[t,var]  = ode45(odefun,tspan,initcond(:,1),opts);
initcond(:,1) = var(end,:);
ind_T = 2;

Var = var(1:end-1,:);
time = t(1:end-1);
Rvec = R(1)*ones(size(time));
figure;
set(gca,'FontSize',20)

subplot(7,13,[60:65,73:78,86:91])    %subplot change here

plot([0 100],[Rstar,Rstar],'k:','LineWidth',1)
hold on
plot([0 100],[Rend,Rend],'k:','LineWidth',1)

if R(1) <= (Rend+Rstar)/2
    plot(...
        [0,T(1)],...
        [R(1),R(1)],...
        'Color',[0 0.6 1],'LineWidth',2);
else
    plot(...
        [0,T(1)],...
        [R(1),R(1)],...
        'Color',[1 0 1],'LineWidth',2);
end
hold on
% plotting the climate
while sum(T)<Tend
    T(ind_T) = nbininv(rand,RR,PnBin);
    while (T(ind_T) == 0)
        T(ind_T) = nbininv(rand,RR,PnBin);
    end
    R(ind_T) = Rend + rand*(Rstar - Rend);
    odefun   = @(t,var)RModefun(var,R(ind_T));
    tspan    = [sum(T(1:ind_T - 1)),sum(T(1:ind_T))];
    [t,var]  = ode45(odefun,tspan,initcond(:,ind_T-1),opts);
    initcond(:,ind_T) = var(end,:);
    
    Var = [Var;var(1:end-1,:)];
    time = [time;t(1:end-1)];
    Rvec = [Rvec;R(ind_T)*ones(size(t(1:end-1)))];
    
    if R(ind_T) <= (Rend+Rstar)/2
        plot(...
            [sum(T(1:ind_T-1)),sum(T(1:ind_T))],...
            [R(ind_T),R(ind_T)],...
            'Color',[0 0.6 1],'LineWidth',3);
    else
        plot(...
            [sum(T(1:ind_T-1)),sum(T(1:ind_T))],...
            [R(ind_T),R(ind_T)],...
            'Color',[1 0 1],'LineWidth',3);
    end
    plot(...
        [sum(T(1:ind_T-1)),sum(T(1:ind_T-1))],...
        [R(ind_T-1),R(ind_T)],...
        'k:','LineWidth',1)
    
    ind_T = ind_T + 1;
end

xticks([0 20 40 60 80 100]);
xticklabels([0 20 40 60 80]);
yticks([1.5 2 2.5 3]);
yticklabels({1.5 '' 2.5});
xlabel('$t$','Position',[98.0358,1.15,-1])
ylabel('$r$','Rotation',0,'Position',[-4.464285512055667,2.849354673773695,-1])
set(gca,'FontSize',20)

annotation('textbox',...
    [0.66875 0.335714285714286 0.1075 0.042857142857143],...
    'Color',[1 0 1],...
    'String','Type-H years',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off');

annotation('textbox',...
    [0.756875 0.137142857142857 0.1075 0.042857142857143],...
    'Color',[0 0.6 1],...
    'String','Type-L years',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off');


axis([0 100 1.2 3.2])

Nvar = @(tt)interp1(time,Var(:,1),tt);
Pvar = @(tt)interp1(time,Var(:,2),tt);
Rvar = @(tt)interp1(time,Rvec,tt);

TPer = [9.665
    9.4381
    9.1755
    11.4994
    10.1756
    10.1489
    8.5541
    8.1481
    11.8877
    12.1880
    9.0007
    10.6828
    10.4563
    10.7188
    8.2345
    8.2154
    9.1937
    11.2690
    10.6334
    9.4462
    8.5651
    11.4770
    11.0905
    8.4176];

%%


% making movie and plotting the phase space
MovieName = 'phase_tipping_movie_RM_nothreshold.mp4';
V = VideoWriter(MovieName,'MPEG-4');
freamsPerSec = 20;
V.FrameRate = freamsPerSec;
open(V);
ind_movie = 1;
tstep = 0.1;
% ttSpan = [0:tstep:55,55*ones(1,10*freamsPerSec),55:tstep:100];
ttSpan = [0:tstep:100];
for tt = ttSpan
    
    % first subplot
    
    subplot(7,13,[60:65,73:78,86:91])    %subplot change here
    
    
    if tt >0
        delete(temp_timeline)
    end
    temp_timeline = plot([tt,tt],[1,3.2],'--k');
    
    % calculate the periodic orbit for R(t)
    indPer = find(abs(Rvar(tt)-R)== min(abs(Rvar(tt)-R)));
    Rper = R(indPer);
    Tper = TPer(indPer);
    initcond = [8   0.01];
    
    ivpfun   = @(t,var)RModefun(var,Rper);
    perfun   = @(t,var,Tper)Tper*RModefun(var,Rper);
    
    [~,tempvar]  = ode45(ivpfun,[0 200],initcond);
    initper  = tempvar(end,:);
    ivpsol = ode45(@(t,var)perfun(t,var,Tper),[0 1],initper, opts);
    
    
    % The boundary condtions for the periodic orbit
    ss = linspace(0,1,200);
    tempinit = @(s)deval(ivpsol,s);
    solinit=bvpinit(ss,tempinit,Tper);
    
    BC=@(WL,WR,Tper)(...
        [WL(1)-WR(1);...
        WL(2)-WR(2);...
        WL(2)-initper(2);... %point phase condition
        ]);
    
    BVPper = bvp5c(perfun,BC,solinit);
    
    
    
    nPoints = 1000;
    ss = linspace(0,1,nPoints);
    
    varper = deval(BVPper,ss);
    
    % the threshold theta
    e2      = [mu,0];
    G       = @(var)RModefun(var,Rper);
    JJ      = MyJacobian(G,e2);
    [eigvic,eigval] = eig(JJ);
    if eigval(2,2)<0
        pert = eigvic(:,2);
    else
        pert = eigvic(:,1);
    end
    pert = 1e-4*pert';
    maninit = e2 + pert;
    [~,varman]   = ode45(ivpfun,[10,0],maninit,opts1);
    
    % second subplot
    
    subplot(7,13,[1:6,14:19,27:32,40:45,53:58,66:71,79:84])    %subplot change here
    
    
    NN = Nvar(0:0.001:tt);
    PP = Pvar(0:0.001:tt);
    
    plot(...
        varper(1,:),mult*varper(2,:),...
        'LineWidth',4,'Color',[0.47 0.67 0.19],'Marker','.')
    hold on
    
    %     plot(...
    %         varman(:,1),mult*varman(:,2),...
    %         'LineWidth',3,'Color','r');
    hold on
    
    
    plot(...
        NN,1000*PP,'-','Color','k')
    
    
    hold on
    if tt< 55
        plot(...
            Nvar(tt),1000*Pvar(tt),'k.',...
            'MarkerSize',30,'LineWidth',1.5)
    else
        plot(...
            Nvar(tt),1000*Pvar(tt),'k.',...
            'MarkerSize',30,'LineWidth',1.5)
    end
    axis([0 15 0 25])
    xticks([0 5 10 15]);
    xticklabels([0 5 10]);
    yticks([0 5 10 15 20 25]);
    yticklabels({'' 5 '' 15 '' ''});
    xlabel('$N$','Position',[14.5,-0.36,-1])
    ylabel('$P$','Rotation',0,'Position',[-0.6,23,-1])
    %     legend(...
    %         '$~~~~\Gamma\big(r(t)\big)$',...
    %         '$~~~~\theta\big(r(t)\big)$',...
    %         '$~\big(N(t),P(t)\big)$',...
    %         'Interpreter','latex','EdgeColor','w')
    legend(...
        '$~~~~\Gamma\big(r(t)\big)$',...
        '$~\big(N(t),P(t)\big)$',...
        'Interpreter','latex','EdgeColor','w')
    hold off
    set(gca,'FontSize',20)
    
    subplot(7,13,[8:13,21:26,34:39])    %subplot change here
    
    plot(...
        0:0.001:tt,NN,'Color',[.6 .6 .6],'LineWidth',3,'Marker','.')
    hold on
    plot(...
        0:0.001:tt,mult*PP,'Color','k','LineWidth',3,'Marker','.')
    axis([0 100 0 25])
    xticks([0 20 40 60 80 100]);
    xticklabels([0 20 40 60 80]);
    yticks([0 5 10 15 20 25]);
    yticklabels({'' 5 '' 15 '' '25'});
    xlabel('$t$','Position',[98.0358,-1.2,-1])
    legend('$~N(t)$','$~P(t)$',...
        'Interpreter','latex','EdgeColor','w')
    hold off
    set(gca,'FontSize',20)
    %
    %     if tt<55
    %         annotation('ellipse',...
    %             [0.379375 0.794285714285714 0.00687500000000002 0.0157142857142857],...
    %             'LineStyle','none',...
    %             'FaceColor','k');
    %     else
    %         annotation('ellipse',...
    %             [0.379375 0.794285714285714 0.00687500000000002 0.0157142857142857],...
    %             'LineStyle','none',...
    %             'FaceColor','k');
    %     end
    
    % if tt<55
    %     annotation('ellipse',...
    %         [0.379375 0.837142857142857 0.00687500000000002 0.0157142857142857],...
    %         'LineStyle','none',...
    %         'FaceColor','k');
    % else
    annotation('ellipse',...
        [0.379375 0.837142857142857 0.00687500000000002 0.0157142857142857],...
        'LineStyle','none',...
        'FaceColor','k');
    %end
    
    
    left   =  .3e3;
    bottom =  0.09e3;
    width  =  1.6e3;
    height =  .7e3;
    
    set(gcf,'Position',[left bottom width height])
    set(gcf,'color','w');
    drawnow
    
    F(ind_movie) = getframe(gcf);
    writeVideo(V,F(ind_movie));
    ind_movie = ind_movie + 1;
end
close(V);
winopen(MovieName)

%% Functions

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

% function [value, isterminal, direction] = myEvent(t,y)
% TOL = 1e-7;
% value      = norm(y)<TOL;
% isterminal = 1;   % Stop the integration
% direction  = 0;   % approch zero from either diractions?
% end