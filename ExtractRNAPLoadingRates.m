function [LoadingRate] = ExtractRNAPLoadingRates(Time,Trace)
% Yang Joon Kim
% This is a function to infer the instantaneous RNAP loading rates from
% MS2 traces (either MeanVectorAP or single traces)
% In future, I need to think about errors, like using bootstrapping, etc.

% Input : MS2 trace(for one cycle)

% Assumptions
% Gene length : hbP2-evePr-MS2.V5-lacZ-alphaTub3'UTR = 
% Elongation rate : 1.54kb/min (Garcia,2013)

Tpeak = 4.5/1.54; %(min)

% Define the time
Time;
dt = median(diff(Time));
TpeakIndex = floor(Tpeak/dt);

% Extract the rate from the MS2 traces
Rate = zeros(size(Time));
for i=2:length(Time)
    if Time(i)<Tpeak
        Rate(i-1) = (Trace(i+1)-Trace(i-1))/(2*dt);
    elseif Time(i)>=Tpeak
        Rate(i-1) = (Trace(i)-Trace(i-1)+Trace(i-TpeakIndex))/dt;
    end
end

LoadingRate = Rate;
end