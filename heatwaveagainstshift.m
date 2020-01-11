function HW=heatwaveagainstshift(Temp,Nyr,nsd, Mn, SD, COLD)
%HEATWAVEAGAINSTSHIFT--heatwave against a shifting baseline
%
% HW=heatwaveagainstshift(Temp,Nyr,nsd {,Mn, SD,COLD})
%
% Temp = m-by-n-by-p array of temperatures.  Temp(:,:,j) contains the
%        temperatures in year j
% Nyr =  number of years to include the rolling reference climatology
% nsd =  number of standard deviations that define a heatwave
% Mn = (optional) m-by-n-by-p array of precomputed means 
% SD = (optional) m-by-n-by-p array of precomputed standard deviations
% COLD=(optional) flag. If COLD evaluates to true, then will find cold
% waves
%
% HW will be a logical array the same size as Temp with "True" if the
% temperature in that year exceeded the prior Nyr's mean by nsd
%

sz=size(Temp);
if(length(sz)==2)
    error('expect Temp to be 3D');
end
HW=nans(size(Temp));

m=sz(3);
I=(1:Nyr)';

HAVEMN=0;
if(nargin>3)
    %have Mn, check if we have SD
    if(nargin==3)
        error('must provide both Mn and SD or none');
    end
    HAVEMN=1;
    if(length(size(Mn))~=3)
        error('Mn must be 3D');
    end
    if(size(Mn,1)~=sz(1) || size(Mn,2)~=sz(2) || size(Mn,3)~=sz(3))
        error('Mn must be the same size as Temp');
    end
    if(length(size(SD))~=3)
        error('SD must be 3D');
    end
    if(size(SD,1)~=sz(1) || size(SD,2)~=sz(2) || size(SD,3)~=sz(3))
        error('SD must be the same size as Temp');
    end
end
if(nargin<6)
    COLD=0;
else
    if(COLD)
        disp('finding cold waves');
    else
        disp('finding heatwaves');
    end
end


for j=1:sz(1);
    for p=1:sz(2);
        dat=squeeze(Temp(j,p,:));
        if(sum(~isnan(dat))>30)
            for k=(Nyr+1):m
                I=(1:Nyr)+(k-Nyr);
                if(~HAVEMN)
                    mn=mean(dat(I));
                    sd=std(dat(I));
                else
                    mn=Mn(j,p,k);
                    sd=SD(j,p,k);
                end
                if(~COLD)
                    HW(j,p,k)=Temp(j,p,k)>mn+nsd*sd;
                else
                    HW(j,p,k)=Temp(j,p,k)<mn-nsd*sd;
                end
            end
        end
    end
end