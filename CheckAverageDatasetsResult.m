function CheckAverageDatasetsResult
% This is a function to check if the AverageDatasets is doing right.

%% 1) Check FractionON average (FractionON from individual embryos)
for i=1:length(FractionON)
    clf
    hold on
    for j=1:length(FractionON(1,1,:))
        plot(0:0.025:1,FractionON(i,:,j))
    end
    pause
end

%% 2) Check MeanVectorAP average
for i=1:length(MeanVectorAP)
    clf
    hold on
    for j=1:length(MeanVectorAP_forPlot(1,1,:))
        plot(0:0.025:1,MeanVectorAP_forPlot(i,:,j))
    end
    pause
end
end