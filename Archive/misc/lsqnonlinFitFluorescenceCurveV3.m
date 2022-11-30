function chi2=lsqnonlinFitFluorescenceCurveV2(TimeData,FluoData,Delay,ncLength,x0)

%V3: Now fits the downward slope independently
% YJK : V3 is for the case of three slopes, 1) Initiation(RNAP
% loading only (Rate0), 2) either plateau or decrease since the loading rate(Rate1) can
% be smaller than the falling off rate(RateOff0), 3) falling off of RNAP only

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
TimeStart0=x0(1);
TimeEnd0=x0(2);
Rate0=x0(3);
RateOff0=x0(4);
Rate1 = x0(5);
% Rate1 is the RNAPII loading rate during the tilted phase.
% RateOff0 is assumed to be negative.

for i=1:length(TimeData)
    if TimeData(i)<=TimeStart0
        FluoPrediction(i)=0;
    elseif (TimeData(i)>TimeStart0)&(TimeData(i)<=TimeStart0+Delay)
        FluoPrediction(i)=Rate0*(TimeData(i)-TimeStart0);
    elseif (TimeData(i)>TimeStart0+Delay)&(TimeData(i)<=TimeEnd0)
        FluoPrediction(i)=(Rate1+RateOff0)*(TimeData(i)-TimeStart0)+Rate0*Delay;
    elseif (TimeData(i)>TimeEnd0)&(TimeData(i)<=(-(Rate0*Delay+(Rate1+RateOff0)*(TimeEnd0-TimeStart0-Delay))/RateOff0+TimeEnd0))
        FluoPrediction(i)=(Rate0*Delay+(Rate1+RateOff0)*(TimeEnd0-TimeStart0-Delay))+...
                            RateOff0*(TimeData(i)-TimeEnd0);
    elseif TimeData(i)>(-(Rate0*Delay+(Rate1+RateOff0)*(TimeEnd0-TimeStart0-Delay))/RateOff0+TimeEnd0)
        FluoPrediction(i)=0;
    end
end



%[TimeRange,Fluorescence]=FluorescenceCurve(ncLength,TimeStart0,TimeEnd0,Rate0,Delay);


%Interpolate the values of the generated curve for the actual time points
%we have
%InterpData=pchip(TimeRange,Fluorescence,TimeData);

%NanFilter=~isnan(FluoData);
%chi2=(FluoData(NanFilter)-InterpData(NanFilter)').^2;
try
    chi2=(FluoData-FluoPrediction').^2;
catch
    1+1;
end

