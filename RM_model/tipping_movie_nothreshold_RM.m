% tipping_movie_RM preduces a movie to illustrate phase-sensetiv tipping in
% the RM model

warning off
clc
clear
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
MovieName = 'phase_tipping_movie_RM_v1_nothreshold.mp4';
V = VideoWriter(MovieName,'MPEG-4');
freamsPerSec = 20;
V.FrameRate = freamsPerSec;
open(V);
ind_movie = 1;
tstep = 0.1;
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
    
    el_1  = 2000;
    el_2  = 2200;
    el_3  = 2400;
    el_4  = 2600;
    el_5  = 2800;
    el_6  = 3000;
    el_7  = 3200;
    el_8  = 3400;
    el_9  = 3600;
    el_10 = 3800;
    
    
    lineWidth  = 1.75;
    markerSize = 3;
    
    linColor1  = [0 0 0];
    linColor2  = [0.1 0.1 0.1];
    linColor3  = [0.2 0.2 0.2];
    linColor4  = [0.3 0.3 0.3];
    linColor5  = [0.4 0.4 0.4];
    linColor6  = [0.5 0.5 0.5];
    linColor7  = [0.6 0.6 0.6];
    linColor8  = [0.7 0.7 0.7];
    linColor9  = [0.8 0.8 0.8];
    linColor10 = [0.9 0.9 0.9];
    
    
    if length(NN)<= el_1
       plot(...
            NN,1000*PP,'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1);
        % #1
    elseif and(length(NN)>el_1,length(NN)<=el_2)
        plot(...
            NN(1:end-el_1),1000*PP(1:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #2
    elseif and(length(NN)>el_2,length(NN)<=el_3)
        plot(...
            NN(1:end-el_2),1000*PP(1:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #3
    elseif and(length(NN)>el_3,length(NN)<=el_4)
        plot(...
            NN(1:end-el_3),1000*PP(1:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #4
    elseif and(length(NN)>el_4,length(NN)<=el_5)
        plot(...
            NN(1:end-el_4),1000*PP(1:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #5
    elseif and(length(NN)>el_5,length(NN)<=el_6)
        plot(...
            NN(1:end-el_5),1000*PP(1:end-el_5),'-','Color',linColor6,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor6,'MarkerEdgeColor',linColor6)
        plot(...
            NN(end-el_5:end-el_4),1000*PP(end-el_5:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #6
    elseif and(length(NN)>el_6,length(NN)<=el_7)
        plot(...
            NN(1:end-el_6),1000*PP(1:end-el_6),'-','Color',linColor7,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor7,'MarkerEdgeColor',linColor7)
        plot(...
            NN(end-el_6:end-el_5),1000*PP(end-el_6:end-el_5),'-','Color',linColor6,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor6,'MarkerEdgeColor',linColor6)
        plot(...
            NN(end-el_5:end-el_4),1000*PP(end-el_5:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #7
    elseif and(length(NN)>el_7,length(NN)<=el_8)
        plot(...
            NN(1:end-el_7),1000*PP(1:end-el_7),'-','Color',linColor8,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor8,'MarkerEdgeColor',linColor8)
        plot(...
            NN(end-el_7:end-el_6),1000*PP(end-el_7:end-el_6),'-','Color',linColor7,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor7,'MarkerEdgeColor',linColor7)
        plot(...
            NN(end-el_6:end-el_5),1000*PP(end-el_6:end-el_5),'-','Color',linColor6,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor6,'MarkerEdgeColor',linColor6)
        plot(...
            NN(end-el_5:end-el_4),1000*PP(end-el_5:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #8
    elseif and(length(NN)>el_8,length(NN)<=el_9)
        plot(...
            NN(1:end-el_8),1000*PP(1:end-el_8),'-','Color',linColor9,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor9,'MarkerEdgeColor',linColor9)
        plot(...
            NN(end-el_8:end-el_7),1000*PP(end-el_8:end-el_7),'-','Color',linColor8,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor8,'MarkerEdgeColor',linColor8)
        plot(...
            NN(end-el_7:end-el_6),1000*PP(end-el_7:end-el_6),'-','Color',linColor7,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor7,'MarkerEdgeColor',linColor7)
        plot(...
            NN(end-el_6:end-el_5),1000*PP(end-el_6:end-el_5),'-','Color',linColor6,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor6,'MarkerEdgeColor',linColor6)
        plot(...
            NN(end-el_5:end-el_4),1000*PP(end-el_5:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #9
    elseif and(length(NN)>el_9,length(NN)<=el_10)
        plot(...
            NN(1:end-el_9),1000*PP(1:end-el_9),'-','Color',linColor10,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor10,'MarkerEdgeColor',linColor10)
        plot(...
            NN(end-el_9:end-el_8),1000*PP(end-el_9:end-el_8),'-','Color',linColor9,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor9,'MarkerEdgeColor',linColor9)
        plot(...
            NN(end-el_8:end-el_7),1000*PP(end-el_8:end-el_7),'-','Color',linColor8,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor8,'MarkerEdgeColor',linColor8)
        plot(...
            NN(end-el_7:end-el_6),1000*PP(end-el_7:end-el_6),'-','Color',linColor7,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor7,'MarkerEdgeColor',linColor7)
        plot(...
            NN(end-el_6:end-el_5),1000*PP(end-el_6:end-el_5),'-','Color',linColor6,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor6,'MarkerEdgeColor',linColor6)
        plot(...
            NN(end-el_5:end-el_4),1000*PP(end-el_5:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
        % #10
    else
        plot(...
            NN(end-el_10:end-el_9),1000*PP(end-el_10:end-el_9),'-','Color',linColor10,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor10,'MarkerEdgeColor',linColor10)
        plot(...
            NN(end-el_9:end-el_8),1000*PP(end-el_9:end-el_8),'-','Color',linColor9,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor9,'MarkerEdgeColor',linColor9)
        plot(...
            NN(end-el_8:end-el_7),1000*PP(end-el_8:end-el_7),'-','Color',linColor8,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor8,'MarkerEdgeColor',linColor8)
        plot(...
            NN(end-el_7:end-el_6),1000*PP(end-el_7:end-el_6),'-','Color',linColor7,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor7,'MarkerEdgeColor',linColor7)
        plot(...
            NN(end-el_6:end-el_5),1000*PP(end-el_6:end-el_5),'-','Color',linColor6,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor6,'MarkerEdgeColor',linColor6)
        plot(...
            NN(end-el_5:end-el_4),1000*PP(end-el_5:end-el_4),'-','Color',linColor5,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor5,'MarkerEdgeColor',linColor5)
        plot(...
            NN(end-el_4:end-el_3),1000*PP(end-el_4:end-el_3),'-','Color',linColor4,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor4,'MarkerEdgeColor',linColor4)
        plot(...
            NN(end-el_3:end-el_2),1000*PP(end-el_3:end-el_2),'-','Color',linColor3,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor3,'MarkerEdgeColor',linColor3)
        plot(...
            NN(end-el_2:end-el_1),1000*PP(end-el_2:end-el_1),'-','Color',linColor2,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor2,'MarkerEdgeColor',linColor2)
        plot(...
            NN(end-el_1:end),1000*PP(end-el_1:end),'-','Color',linColor1,'LineWidth',lineWidth,...
            'Marker','.','MarkerSize',markerSize,'MarkerFaceColor',linColor1,'MarkerEdgeColor',linColor1)
    end
    hold on
    
    plot(...
        Nvar(tt),1000*Pvar(tt),'k.',...
        'MarkerSize',30,'LineWidth',1.5)
    
    axis([0 15 0 25])
    xticks([0 5 10 15]);
    xticklabels([0 5 10]);
    yticks([0 5 10 15 20 25]);
    yticklabels({'' 5 '' 15 '' ''});
    xlabel('$N$','Position',[14.5,-0.36,-1])
    ylabel('$P$','Rotation',0,'Position',[-0.6,23,-1])
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
    
annotation('ellipse',...
    [0.379375 0.837142857142857 0.00687500000000002 0.0157142857142857],...
    'LineStyle','none',...
    'FaceColor',[0 0 0]);

annotation('line',[0.359975 0.379975],[0.844328571428571 0.844328571428571],'LineWidth',2);
    
    
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

function df = MyJacobian(f,x)
% to coumpute a jacobian matrix for a function f:R^n --> R^n
h=1e-4;
n=length(f(x)); m=length(x);
df = NaN(n,m);
if isfinite(sum(x))
    % F=repmat(f(x),1,m);
    for j=1:m
        F1=f(x);
        x(j)=x(j)+h;
        F2=f(x);
        x(j)=x(j)-h;
        df(:,j)=(F2-F1)./h;
    end
else
    df = NaN(n,m);
end
end% tipping_movie_RM preduces a movie to illustrate phase-sensetiv tipping in
% the RM model

