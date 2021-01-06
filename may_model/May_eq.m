function [Neq] = May_eq(R,q)
%May_eq is a function of $R$ to find the 3 nontriveal equlibria of May
%model 

C       = 1.75/8;    %/yr.               Nonlinear death rate of prey
alpha   = 505;       % prey/(pred.yr)    Predator saturating kill rate
beta    = 0.3;       % prey/ha           Predator half-saturating constant
s       = 0.85;      % /yr                  Rate of predator population.
mu      = 0.03;      % NA                 Allee parameter.
nu      = 0.003;      % NA                Allee parameter.
eps     = 0.031;       % NA                Artificial parameter.


% poly = [-q*C,...
%     q*R - q*C*beta + q*C*mu - alpha,...
%     q*R*beta + q*C*mu*beta - q*R*mu - alpha*eps - alpha*nu,...
%     -q*R*beta*mu - alpha*eps*nu];

poly1 = [C*q,...
    -(q.*R + C.*q.*mu - beta.*C.*q -alpha ),...
    -(beta.*q.*R + beta.*C.*q.*mu - q.*R.*mu - alpha.*(eps + nu)),...
    (beta.*q.*R.*mu + alpha.*nu.*eps)];

Ntemp = roots(poly1);

if and(imag(Ntemp(1))==0,real(Ntemp(1))>=0)
    Neq(1) = Ntemp(1);
else
    Neq(1) = NaN;
end

if and(imag(Ntemp(2))==0,real(Ntemp(2))>=0)
    Neq(2) = Ntemp(2);
else
    Neq(2) = NaN;
end

if and(imag(Ntemp(3))==0,real(Ntemp(3))>=0)
    Neq(3) = Ntemp(3);
else
    Neq(3) = NaN;
end
% Neq = Ntemp;
end

% test 
% P =(N+eps)./q;
% R.*N.*(1-(C/R).*N).*((N-mu)./(nu+N)) - alpha.*N.*P./(N+beta)


%% the drevation 

% R     = sym('R'); 
% C     = sym('C');
% alpha = sym('alpha');    
% beta  = sym('beta');
% s     = sym('s');     
% q     = sym('q'); 
% mu    = sym('mu');  
% nu    = sym('nu'); 
% eps   = sym('epslion');
% 
% N = sym('N');
% 
% result1 = ( ( q.*R-((C.*q).*N) ).*(N-mu).*(N+beta) ) + ...
%     (-alpha.*(N+eps).*(nu+N));
% 
% result2 = ( q.*R.*N - q.*R.*mu - C.*q.*(N.^2) + C.*q.*mu.*N ).*...
%     (N+beta) + (-alpha.*(N+eps).*(nu+N));
% 
% result3 = ( (q.*R.*(N.^2) - q.*R.*mu.*N - C.*q.*(N.^3) + C.*q.*mu.*(N.^2) +...
%      beta.*q.*R.*N - beta.*q.*R.*mu - beta.*C.*q.*(N.^2) + beta.*C.*q.*mu.*N )) +...
%     (-alpha.*(N.^2) + N.*(-alpha*eps-alpha*nu) - alpha*nu*eps);
% 
% result4 = (N.^3) .* (- C.*q) + ...
%     (N.^2).*(q.*R + C.*q.*mu - beta.*C.*q -alpha ) +...
%     N.*(beta.*q.*R + beta.*C.*q.*mu - q.*R.*mu -alpha.*eps - alpha.*nu) + ...
%     (-beta.*q.*R.*mu - alpha*nu*eps);
% 
% poly = [-q*C,...
%     q*R - q*C*beta + q*C*mu - alpha,...
%     q*R*beta + q*C*mu*beta - q*R*mu - alpha*eps - alpha*nu,...
%     -q*R*beta*mu - alpha*eps*nu];
% 
% Eq = roots(poly);
%%
% C     = 1.75/8;    %/yr.               Nonlinear death rate of prey
% alpha = 505;       % prey/(pred.yr)    Predator saturating kill rate
% beta  = 0.3;       % prey/ha           Predator half-saturating constant
% s     = 0.85;      % /yr               Rate of predator population.
% q     = 205;       % prey/pred         Minimum prey biomass.
% mu    = 0.03;      % NA                Allee parameter.
% nu    = 0.003;      % NA                Allee parameter.
% eps   = 0.031;
% pp = @(R,N)(-q*C*N^3 + ...
%     (q*R - q*C*beta + q*C*mu - alpha)*N^2 + ...
%     (q*R*beta + q*C*mu*beta - q*R*mu - alpha*eps - alpha*nu)*N + ...
%     -q*R*beta*mu - alpha*eps*nu);
%% The model with out Allee 

% global  C alpha beta s q mu nu eps
% 
% R     = 4 %sym('R');     %/yr.            Prey intrinsic growth rate
% C     = 1.75/8;    %/yr.               Nonlinear death rate of prey
% alpha = 505;       % prey/(pred.yr)    Predator saturating kill rate
% beta  = 0.3;       % prey/ha           Predator half-saturating constant
% s     = 0.85;      % /yr               Rate of predator population.
% q     = 212;       % prey/pred         Minimum prey biomass.
% mu    = 0.2;      % NA                Allee parameter.
% nu    = 0.72;      % NA                Allee parameter.
% eps   = 0.1;
% 
% 
% poly = [-q*C,...
%     q*R-q*C*beta- alpha,...
%     q*R*beta-alpha*eps];
% 
% N=roots(poly);
% % q*R.*Eq.*(1-(C/R).*Eq)-(alpha.*(Eq.^2))./(beta+Eq)
% q.*R.*(1-(C/R).*N).*(beta+N) - alpha.*(N+eps)