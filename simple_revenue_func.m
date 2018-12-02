function R=simple_revenue_func(t, tmn, sd)
%SIMPLE_REVENUE_FUNC--simple function that computes $$ given temperature
%
% R=simple_revenue_func(t, tmn, sd)
%
% t = temperature (or other independent variable)
% tmn = temperature where revenue is maxed
% sd = standard deviation of temperature
%
% R = normpdf(t,tmn,sd)-normpdf(tmn+2*sd,tmn,sd)
%
%  R is >=0 if abs(t-tmn)<2*sd
%


R=1/sqrt(2*pi*sd*sd)*(exp(-(t-tmn).^2/(2*sd*sd))-exp(-2));
