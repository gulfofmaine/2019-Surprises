function varargout=evaluate_strategy_over_time(params,algorithm,M,envinfo,switchcost, discount)
%EVALUATE_STRATEGY_OVER_TIME--evaluate a strategy in an environment
%
%     Net = evaluate_strategy_over_time(params,algorithm,M,envinfo,switchcost,discount)
% or
%   [Net, RCN]=evaluate_strategy(...)
% or
%   [Net, RCN, tBT]=evaluate_strategy(...)
% 
% params = [w, a, b]--weighting for the different choices.
%         w = weight placed on not losing
%         a = weight placed on prior revenue when deciding whether to
%             switch
%         b = resistance to change
% algorithm = 'back' or 'trend'--algortihm for figuring out the environment
% M = time window (years)
% envinfo = either m-by-1 list of temperatures (m>M) or
%           m-by-3 list of environmental parameters where
%           envinfo(j,:)=[rate, gamma, nyrs].  In this case, the
%           mean temperature will be increased by rate each year and the
%           realized temperature will be drawn from a distribution with SD
%           gamma.  nyrs is the number of years to apply these stats before
%           switching to the next stats.  Note that envinfo(1,:) is applied
%           for the M spinup years and then for its nyrs.
% switchcost = cost of switching
% discount = discount rate
%
% The environment will spin up for M years.  Then, we compute the
% statistics (mean or trend), then find the t0 and beta that maximize the
% value (expected return + w*Prb{no loss}). This will be the initial
% strategy. Temperature will be computed (or read if envinfo contains
% actual temperatures, revenue computed.  After than, in each year, the
% statistics will be updated and we will decide whether to switch
% strategies. The algorithm will be run for q=sum(envinfo(:,3)) years if
% temperature is computed or q= m-M years if provided directly
%
% Net = net income over the entire simulation (discounted)
% RCN = q-by-3 array of [revenue, cost, net] in each simulated year
% tBT = q+M-by-3 array of [assumed mean, beta, realized temperature].  Note
%            that tBT(1:M,1:2) = nan
%
% Andrew Pershing (apershing@gmri.org), 2018

if(length(params)~=3)
    error('params must be length 3 = [w, a, b]');
end


algorithmnames={'back','trend'};
pstrat=validatestring(algorithm,algorithmnames);%pstrat is one of the pfunc names

[m,n]=size(envinfo);
REALDATA=0;
if(n==1)
    REALDATA=1;
    if(m<M)
        error('not enough data. envinfo has %d temperatures but window is M=%d\n',m,M);
    end
    q=m-M;
elseif(n==3)
    q=sum(envinfo(:,3));%number of years
else
    error('envinfo should be m-by-1 (real temperatures) or m-by-3 (computed temperatures');
end

%0. initilize RCN

RCN = zeros(q,3);

%1. fill in the temperatures for the spinup period (1:M) and the simulation
%   period (M+1:M+q)
tBT=nans(q+M,3);
if(REALDATA)
    tBT(:,3)=envinfo;
    kst=M+1;
else
    tBT(:,3)=randn(q+M,1);
    tBT(1:M,3)=tBT(1:M,3)*envinfo(1,2)+(0:M-1)'*envinfo(1,1);%computed temperatures
    jj=M+1;
    T0=M*envinfo(1,1);
    for k=1:m
        for j=1:envinfo(k,3);%years to apply this strategy
            T0=T0+envinfo(k,1);
            tBT(jj,3)=tBT(jj,3)*envinfo(k,2)+T0;
            jj=jj+1;
        end
    end
end

%2. get expected temperature and SD from the algorithm
[T0,Gamma, R]=feval(pstrat,tBT(1:M,3));

%3. get our starting strategy
j=M+1;
tBT(j,1)=T0;
      %optimal_strategy_with_discounting(Test, r, gamma,w,M, discount)
TBopt=optimal_strategy_with_discounting(T0,R, Gamma,params(1),M,discount);%optimal beta allowing for discounting
tBT(j,1:2)=TBopt(:)';
jj=1;

npv=1;
for jj=1:q
    %get the realized revenue based on current temperature and strategy
    RCN(jj,1)=simple_revenue_func(tBT(j,3),tBT(j,1),tBT(j,2));
    
    if(jj<q)
        %update the statistics
        [T0,Gamma, R]=feval(pstrat,tBT(j-M+1:j,3));

        %potential new strategy
        TBopt=optimal_strategy_with_discounting(T0,R, Gamma,params(1),M,discount);%optimal beta allowing for discounting

        if(abs(TBopt(1)-tBT(j,1))>1e-6);%T0 sufficiently different
            Ttest=linspace(tBT(j,1),TBopt(1),5);
        else
            Ttest=tBT(j,1);
        end

        if(abs(TBopt(2)-tBT(j,2))>1e-6);%Bopt sufficiently different
            Btest=linspace(tBT(j,2),TBopt(2),5);
            Btest=Btest(:);
        else
            Btest=tBT(j,2);
        end

        mm=length(Btest);
        nn=length(Ttest);

        if(nn>1 || mm>1)
            Vcand=zeros(mm,nn);
            Ccand=Vcand;
            %Cost = integral((N(t|tBT(j,1),tBT(j,2))-N(t|Ttest,Btest)^2 dt)
            %     =integral((N1-N2)^2)
            %     =integral(N1*N1)-2integral(N1*N2)+integral(N2*N2)
            Ccand(:)=integral_N1_N2(tBT(j,1),tBT(j,2),tBT(j,1),tBT(j,2));%integral(N1*N1)

            for rr=1:mm
                for cc=1:nn
                               % expected_revenue_over_time2(t0,beta,Test, r,gamma,w, M,discount)
                    Vcand(rr,cc)=expected_revenue_over_time2(Ttest(cc), Btest(rr), T0, R, Gamma, params(1),M,discount);

                    Ccand(rr,cc)=Ccand(rr,cc)+integral_N1_N2(Ttest(cc), Btest(rr),Ttest(cc), Btest(rr))+...%integral(N2*N2))
                                 -2*integral_N1_N2(Ttest(cc), Btest(rr),tBT(j,1),tBT(j,2));%integral(N1*N2))
                end
            end
            Ccand=Ccand*switchcost;
            
            Vcand=Vcand-Vcand(1,1);%change in value

            [Vbest,I]=max(Vcand(:)-Ccand(:));
            if(Vbest>0 & Vbest>(params(2)*RCN(jj,1)+params(3)))
                %switch
                [rr,cc]=ind2sub([mm,nn],I(1));
                RCN(jj,2)=Ccand(rr,cc);
                tBT(j+1,1:2)=[Ttest(cc),Btest(rr)];
            else
                tBT(j+1,1:2)=tBT(j,1:2);
            end
        else
            tBT(j+1,1:2)=tBT(j,1:2);%same strategy
        end
    end
    
    j=j+1;
end

RCN(:,3)=RCN(:,1)-RCN(:,2);

Net=sum(RCN(:,3));%sum of the net revenue

if(nargout<=1)
    varargout={Net};
elseif(nargout==2)
    varargout={Net,RCN};
elseif(nargout>=3)
    varargout={Net,RCN,tBT};
end



function [T0,Gamma, R]=back(T)
T0=mean(T);
Gamma=std(T);
R=0;



function [T0,Gamma,R]=trend(T)

N=length(T);
Y=(0:N-1)';
coefs=[Y,ones(N,1)]\T;%compute the trend
r=coefs(1);
sigma=std(T);
mn=mean(T);
%Gamma = sigma.^2+(r.^2)/12*(N^2-10); %variance--old forumula
Gamma = (2*N)/(2*N-1)*(sigma.^2-(1/12)*(N^2+N-7)*r^2); %new forumula
if(Gamma<0)
    error('Problem estimating gamma: N=%d, r=%f, sigma=%f',N,r,sigma);
end
Gamma = sqrt(Gamma);
T0=mn+r*(N-(N-1)/2);
R=r;





