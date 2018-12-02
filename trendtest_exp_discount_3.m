%trendtest_exp_discount_3.m
%
% Runs the economic model for the paper
%
% Andrew Pershing (apershing@gmri.org), 2018
%
w=[0,1,2];%weight placed on not losing money
switchcost=[0 1 10 100];%cost of switching
r=linspace(0,0.1,20);%warming trend
gamma=0.1:0.1:1;%temperature variance
Bknet=nans(length(r),length(gamma),length(w),length(switchcost));
Trnet=Bknet;

discount=0.03;%discount rate

for p=1:1;
    for q=length(switchcost):-1:1;fprintf('%2d,%2d\n',p,q);
        for j=1:length(r);fprintf('\t%2d\n',j);
            for k=1:length(gamma);
                envinfo=[r(j),gamma(k),100];
                for i=1:10;
                    [Net,RCNbk,tBTbk]=evaluate_strategy_over_time([w(p),0,0],'back',30,envinfo,switchcost(q),discount);
                    [Net,RCNtr,tBTtr]=evaluate_strategy_over_time([w(p),0,0],'trend',30,tBTbk(:,3),switchcost(q),discount);
                    NetR(i,:)=[sum(RCNbk(:,3)),sum(RCNtr(:,3))];
                    CostR(i,:)=[sum(RCNbk(:,2)),sum(RCNtr(:,2))];
                    if(isnan(NetR(i,1)) || isnan(NetR(i,2)))
                        disp('Oh heck. Nans');
                    end
                end
                Bknet(j,k,p,q)=mean(NetR(:,1));Bkcost(j,k,p,q)=mean(CostR(:,1));
                Trnet(j,k,p,q)=mean(NetR(:,2));Trcost(j,k,p,q)=mean(CostR(:,2));
            end;
        end;
        save TrendTest_discounted_running4 switchcost w r gamma *net *cost discount
    end;
end