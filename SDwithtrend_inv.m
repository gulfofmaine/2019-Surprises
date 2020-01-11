function gamma=SDwithtrend_inv(r,N,sigma)
%SDWITHTREND_INV--compute std. deviation if there is a trend
%
% gamma=SDwithtrend_inv(r,N,sigma)
%
% The variance assuming stationarity (sigma^2) involves both the interannual 
% variance (gamma^2) and the trend (r). The weighting of these two terms
% depends on the number of years (N).
%
% This function takes the total variance (sigma) and computes the
% interannual varicne (gamma).
%
% 
coefsf1=[-1/2
1];
A=coefsf1(1)/N+coefsf1(2);

coefsN=[1/12
1/12
-7/12];
B=[N^2,N,1]*coefsN;

gamma=(sigma.^2-B*r.^2)/A;
gamma=sqrt(gamma);