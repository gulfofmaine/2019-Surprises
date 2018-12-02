function varargout=expected_revenue_over_time2(t0,beta,Test, r,gamma,w, M,discount)
%EXPECTED_REVENUE_OVER_TIME2
%
% D=expected_revenue_over_time2(t0,beta,r,Test,gamma,w, M,discount)
%
% computes the expected revenue (assuming the simple revenue function) for
% a given investment strategy and in an assumed environment. Revenue is 
% discounted over M years, and the envrionmental temperature is allowed to
% have a linear trend.
%
%   t0 = mean temperature of the investment strategy
% beta = standard deviation of the investment strategy
% Test = temperature in the environment at the beginning
%    r = warming rate
% gamma= standard deviation of the environmental temperature.
%    w = weighting placed on avoiding losses
%    M = time horizon for investments
% discount = discount rate.
%
% WR = R + w*P
%
% R and P are computed using the analytic solution to the integral
%
%
%

j=(0:M-1)';
T=Test+r*j;%temperatures

%R = integral(N(x|T,gamma)*R(x|0,beta),dx)
%  = integral(N(x|T,gamma)*N(x|0,beta),dx)-exp(-2)./sqrt(2*pi*beta*beta)
R=zeros(M,1);
for k=1:M
    R(k)=expected_revenue2(t0, beta,T(k), gamma,w);%net value
end
D=(1+discount).^-j;

tDR=sum(D.*R);%total revenue

switch nargout
    case 1
        varargout={tDR};
    case 2
        varargout={tDR,[D,R]};
    otherwise
        varargout={tDR};
end
