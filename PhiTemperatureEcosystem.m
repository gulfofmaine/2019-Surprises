function [N,T] = PhiTemperatureEcosystem(tau,sigma,r0, mu, phi, nyrs, envinfo, GAMMA,  N0)
%PHITEMPERATUREECOSYSTEM--neutral ecosystem model under a warming trend 
%
% [N, T] = PhiTemperatureEcosystem(tau,sigma,r0, mu, phi, nyrs, envinfo, GAMMA, N0)
%
% tau   = length m vector of temperatures preferences for the m species
% sigma = length m vector of temerpature variances for the m species
% r0    = average population growth rate
% mu     = constant density dependent mortality rate for all species
% phi = the phi parameter from the diversity model
%         0 => competitive exclusion
%         1 => non-interactive coexistence
% envinfo = nyrs+1-by-1 array of temperatures or a struct with fields
%       TEM   = starting temperature in the environment
%       GAM   = interannual variance in temperature
%       SLOPE = warming rate (degrees per year)
% GAMMA = interannual variance in temperature (should be the same as
%       envinfo.GAM if provided)
% N0    = initial abundance.  length m vector (specified abundance for each
%         species), scalar (same abundance for each), or if empty or 
%         missing, then K is useed.
%
% N will be m-by-nyrs+1 and will contain the abundance in each year
% T will be nyrs+1-by-1 and will contain the temperatures
%
% Andrew Pershing (apershing@gmri.org), 2018

if length(tau)~=length(sigma)
    error('tau and sigma should be the same length')
end
tau=tau(:);
sigma=sigma(:);

m=length(tau);

Nstart=zeros(m,1);;

if ~exist('N0')
    Nstart(:) = 1/m; % Initial population size
else
    Nstart(:) = N0(:);
end

%get the temperature in each year
if(isstruct(envinfo))
    T=randn(nyrs+1,1)*envinfo.GAM+envinfo.TEM+envinfo.SLOPE*(0:nyrs)';
    if(envinfo.GAM~=GAMMA)
        error('envionfo.GAM and GAMMA should be the same');
    end
else
    if(length(envinfo)~=nyrs+1)
        error('envinfo should be a struct or have nyrs+1 entries');
    end
    T=envinfo(:);
end

%precompute the coefficients on growth rate to"neutralize" the model
Rcoef=zeros(m,1);
for j=1:m;
    Rcoef(j)=r0/integral_N1_N2(0,GAMMA,0,sigma(j));
end

%set up the ode function
odefunc=@(t,N)phitempODE(t,N,Rcoef,tau,sigma,T, mu,phi);

%make sure we're dealing with positive numbers=
odeopt=odeset('NonNegative',ones(m,1));
[t,N]=ode45(odefunc,0:nyrs,Nstart,odeopt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dN=phitempODE(t,N,Rcoef,tau,sigma,T,mu,phi);
%PHITEMPODE--ODE function with the phi model
%

t0=interp1((0:length(T)-1)',T,t);%linearly interpolate temperature
R=Rcoef;
I=find(N<0);
N(I)=0;
for j=1:length(Rcoef);
    R(j)=Rcoef(j).*normpdf(t0,tau(j),sigma(j));%the growth rate for each species
end
Nbar=sum(N);
dN=R.*N-mu*(N.^(1+phi).*(Nbar.^(1-phi)));
dN(I)=-1e-12*N(I);%a small, positive change to try to nudge back to zero

% I=find(abs(imag(dN)));
% if(~isempty(I));
%     disp('imaginary');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dN=phitempODE_logistic(t,N,Rcoef,tau,sigma,T,Kall,r0, src,phi);
%PHITEMPODE--ODE function with the phi model
%

t0=interp1((0:length(T)-1)',T,t);%linearly interpolate temperature
R=Rcoef;
K=R;
for j=1:length(Rcoef);
    R(j)=Rcoef(j).*normpdf(t0,tau(j),sigma(j));%the growth rate for each species
    K(j)=Kall*R(j)/r0;
end
Nbar=sum(N);
dN=R.*N.*(K-(N.^phi).*(Nbar.^(1-phi)))./K+src;




