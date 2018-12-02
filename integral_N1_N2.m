function S=integral_N1_N2(mn1,sd1,mn2,sd2)
%INTEGRAL_N1_N2 integral of the product of two normal distributions
%
% S=integral_N1_N2(mn1,sd1,mn2,sd2)
%
% mn1, mn2 = means of the two distributions
% sd1, sd2 = SDs of the two distributions
%

C=1/sqrt(2*pi*sd1^2*(1+(sd2/sd1)^2));

X=-(mn1-mn2)^2/(2*sd1^2*(1+(sd2/sd1)^2));

S=C*exp(X);
