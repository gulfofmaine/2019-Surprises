%EcosystemExp_2b.m
%
% Runs the community model experiment reported in the paper
%
% Andrew Pershing (apershing@gmri.org), 2018

%set up the environment
trend=linspace(0,0.1,20);%warming trends
gamma=0.1:0.1:1;%temperature variance

%setup biology
taumax=max(trend)*100+2*max(gamma);%max tau we need to cover the range of environments
tau=linspace(-2*max(gamma),taumax,20);%preferred temperatures
sigma=linspace(0.1,0.75,15)';%temperature variances

tau=repmat(tau(:)',length(sigma),1);%create combinations
sigma=repmat(sigma(:),1,length(tau));

t2=[0.5 1 2 4];%doubling time
r0=log(2)./t2;%growth rates
mu=1e-2;%mortality rates

phi=[0.25 0.5 0.75 1];

reps=10;
for a1=3:4;%phi
    fprintf('%2d/%2d\n',a1,length(phi));
    for a2=1:length(r0);
        
        fprintf('\t%2d/%2d\n',a2,length(r0));
        
        Eco(a1,a2).phi=phi(a1);
        Eco(a1,a2).r0=r0(a2);
        Eco(a1,a2).mu=mu;
        
        Eco(a1,a2).N0=nans(reps,size(tau,1)*size(tau,2),length(trend),length(gamma));
        Eco(a1,a2).NF=nans(reps,size(tau,1)*size(tau,2),length(trend),length(gamma));
        
        for b1=1:length(trend);
            for b2=1:length(gamma);
                for j=1:reps;
                    envinfo=struct('TEM',0,'GAM',gamma(b2),'SLOPE',0);
                    [N0, T]=PhiTemperatureEcosystem(tau(:), sigma(:),r0(a2),mu,phi(a1),500,envinfo,envinfo.GAM,1e-3);
                    envinfo.SLOPE=trend(b1);
                    [N1, T]=PhiTemperatureEcosystem(tau(:), sigma(:),r0(a2),mu,phi(a1),100,envinfo,envinfo.GAM,N0(end,:));
                    Eco(a1,a2).N0(j,:,b1,b2)=N1(1,:);
                    Eco(a1,a2).NF(j,:,b1,b2)=N1(end,:);
                end
            end
        end
        save EcosystemExp_running2b trend gamma tau sigma r0 mu phi Eco
    end
end