function [TotalmRNA_averaged, TotalmRNA_SEM] = Average_TotalmRNAProd(DataType, DataPath)
% DESCRIPTION
% Let's use AccumulatedmRNA_FractionON (from AverageDatasets.m)
% stored for each embryo calculated by AverageDatasets.m

D = dir(DataPath);

k=1; % index for the Datasets
for i=1:length(D)
    if contains(D(i).name,DataType,'IgnoreCase',true)
        DataName{k} = D(i).name;
        k=k+1;
    end
end

% Load processed AccumulatedmRNA_FractionON from individual datasets
% (individual embryos)
for i=1:length(DataName)
    Data = load([DataPath,filesep,DataName{i}])
    AccumulatedmRNA_compiled{i} = Data.AccumulatedmRNA_FractionON;
    % Total mRNA is defined using the value at the last point 
    % (from NC12 to almost end of NC14)
    TotalmRNAProd(:,i) = Data.AccumulatedmRNA_FractionON(end,:);
end

% Check by plotting
% hold on
% plot(0:0.025:1,TotalmRNAProd(:,1))
% plot(0:0.025:1,TotalmRNAProd(:,2))
% plot(0:0.025:1,TotalmRNAProd(:,3))


% Convert zeros to Nans, actually this should be done in the
% AverageDatasets.m
TotalmRNAProd(TotalmRNAProd==0) = nan;

TotalmRNA_averaged = nanmean(TotalmRNAProd,2)
TotalmRNA_SEM = nanstd(TotalmRNAProd,[],2)./sqrt(length(DataName));

%% Plot TotalmRNA (MS2) - Mean and individual - r0
hold on
errorbar(0:0.025:1, TotalmRNA_averaged, TotalmRNA_SEM)
for i=1:length(DataName)
    h(i) = plot(0:0.025:1,TotalmRNAProd(:,i))
end

title(['Accumulated mRNA over AP for ',DataType(1:end-1)])
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('Mean')

StandardFigure(gcf,gca)
% Save this figure
saveas(gcf,[DataPath,filesep,DataType,'AccumulatedmRNA_Mean_Individual.tif'])
saveas(gcf,[DataPath,filesep,DataType,'AccumulatedmRNA_Mean_Individual.pdf'])
clf
end