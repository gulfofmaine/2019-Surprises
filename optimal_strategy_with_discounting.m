function TBopt=optimal_strategy_with_discounting(Test, r, gamma,w,M, discount)
%OPTIMAL_STRATEGY_WITH_DISCOUNTING--optimal beta value for simple_revenue
%
% TBopt=optimal_strategy_with_discounting(Test, r, gamma,w,M, discount)
%
% Finds the strategy [t0, beta] that maximizes the discounted returns in
% the assumed environment.
%
%  Test = starting temperature in the environment
%     r = warming trend
% gamma = standard deviation of the environment
%     w = weighting placed on avoiding losses
%     M = number of years for investment horizon
% discount = discount rate
%
% TBopt will be the [t0, beta] that maximizes 
%
%       WR(beta) = ER(beta) + w* Pnz(beta)
%
% where ER(beta) is the expected revenue and Pnz is the probability of not
% going bust
%
% Andrew Pershing (apershing@gmri.org), 2018




if(r~=0)
             % expected_revenue_over_time2(t0,beta,Test, r,gamma,w, M,discount)
    objf=@(x)(-expected_revenue_over_time2(x(1),x(2),Test, r,gamma,w,M,discount));

    optopt=optimset('fminsearch');
    optopt=optimset(optopt,'Display','off');
    [TBopt,negfval,flag]=fminsearch(objf,[r,gamma],optopt);
    if(flag<0)
        %try one more time
        fprintf('fminsearch flag = %d. Trying again\n',flag);
        I=find(isnan(TBopt)|isinf(TBopt));
        if(length(I)==2)
            TBopt=[1,1];%try something else
        elseif(length(I)==1)
            tmp=[0,0.5];
            TBopt(I)=tmp;
        end
        TBopt=TBopt.*(1+randn(size(TBopt))*1e-2);%add some noise
        [TBopt,negfval,flag]=fminsearch(objf,TBopt,optopt);
    elseif(flag==0)
        [TBopt,negfval,flag]=fminsearch(objf,TBopt,optopt);%try one more time
    end
else
    %constant environment.  No need to discount since
    %the temperature is assumed to never change. Discounting would change
    %the return, but it doesn't affect the optimal strategy
    %
             % expected_revenue2(t0, beta,T, gamma,w)
    objf=@(B)(-expected_revenue2(Test,B,Test,gamma,w));
    Bopt = fminbnd(objf,0.001*gamma,4*gamma);
    TBopt=[Test,Bopt];
end




