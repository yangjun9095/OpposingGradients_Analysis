function ConvertMS2ToRNAPLoadingRate
% This script is for extracting the RNAP loading rate from the MS2 traces

%% Load the dataset
Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3_FromNC12.mat')

%% Define the inputs for ExtractRNAPLoadingRates.m
MeanVector13 = Data.MeanVectorAP(nc13:nc14,:);
Time = Data.ElapsedTime(nc13:nc14)-Data.ElapsedTime(nc13);

%% Extract the Loading Rates for all AP bins
LoadingRate = nan(length(Time),41)
for i=11:26 % for AP bins
    LoadingRate(:,i) = ExtractRNAPLoadingRates(Time,MeanVector13(:,i));
end


%% Plot
hold on
for i=11:26
    plot(Time,LoadingRate(:,i))
    pause
end
end