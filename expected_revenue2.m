function WR=expected_revenue2(t0, beta,T, gamma,w)
%EXPECTED_REVENUE--compute the expected revenue, including loss aversion
%
% WR=expected_revenue2(t0, beta,T, gamma,w)
%
% computes the expected revenue (assuming the simple revenue function) for
% a given investment strategy and in an assumed environment. 
%
%   t0 = mean temperauter of the investment strategy
% beta = standard deviation of the investment strategy
%    T = mean temperature in the environment 
% gamma= standard deviation of the environmental temperature.
%    w = weighting placed on avoiding losses
%
% WR = R + w*P
%
% R and P are computed using the analytic solution to the integral
%

R=integral_N1_N2(T,gamma,t0,beta);%analytic solution of integral
ER=R-exp(-2)./sqrt(2*pi*beta*beta);%

if(w>0)
    %integral of temperature likelihood over the range of temperatures where the investment is >0
    WR=normcdf(t0+2*beta,T,gamma)-normcdf(t0-2*beta,T,gamma);
else
    WR=0;
end
WR=w*WR+ER;