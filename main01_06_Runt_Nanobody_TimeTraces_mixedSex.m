function main01_06_Runt_Nanobody_TimeTraces_mixedSex
% This script is for generating plots of Runt-Nanobody (BG-subtracted)
% time-traces.
% We will use the datasets that are already background-subtracted.

%% Step1. Average the BG-subtracted datasets for synchronization
AverageDatasets_LlamaTaggedProtein('Runt-1min-200Hz-mixed','savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData')

%% Step2. Load the mat file for plotting

FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
Runt_mixed = load([FilePath, filesep, 'Runt-1min-200Hz-mixed_BGsubtracted-Averaged.mat']);

%% Extract useful fields
% Time, MeanVectorAP_indiviual, SDVector_individual,
% NParticlesAP_individual, etc.

% mixed
MeanFluo = Runt_mixed.MeanVectorAP_BGsubtracted_individual;
SDFluo= Runt_mixed.SDVectorAP_BGsubtracted_individual;
NNuclei = Runt_mixed.NParticlesAP_BGsubtracted_individual;
SEFluo = SDFluo./sqrt(NNuclei); % for SEM of individual embryos.
Time = Runt_mixed.ElapsedTime;
NC13 = Runt_mixed.nc13;
NC14 = Runt_mixed.nc14;

%% Plot for checking (over AP axis)
% APaxis = 0:0.025:1;
% EmbryoIndex = 3;
% 
% hold on
% for i=1:length(MeanFluo(:,1,EmbryoIndex))
%     errorbar(0:0.025:1, MeanFluo(i,:,EmbryoIndex), SEFluo(i,:,EmbryoIndex))
%     pause
% end


%% Plot for checking (individual traces over Time) : Check Synchronization
tLength = length(Time);
APbin = 15;
hold on
% female
for embryo=1:length(MeanFluo(1,1,:))
    errorbar(Time(1:tLength), MeanFluo(1:tLength,APbin,embryo), SDFluo(1:tLength,APbin,embryo),'r')
    pause
end

%% Averaging over embryos
% Female
numEmbryos = length(MeanFluo(1,1,:));
Averaged_Fluo = nanmean(MeanFluo, 3);
SD_Fluo = nanstd(MeanFluo, 0 , 3);
SE_Fluo = SD_Fluo./sqrt(numEmbryos); % # of embryos
Time = Time(1:tLength); % New time

%% Plot - Averaged trace over Time
for APbin = 9:25 % Covers 20%-60% of the embryo
    clf
    APpos = (APbin - 1) * 2.5;

    hold on
    errorbar(Time, Averaged_Fluo(1:tLength,APbin), SE_Fluo(1:tLength,APbin))

    % Vertical bar for NC demonstration
    xline(NC14,'--k')
    % xlim([0 max(Time_female)])
    xlim([0 50])
    ylim([100 max(Averaged_Fluo(1:50,APbin)) + 50])
    title('Runt protein concentration over Time')
    xlabel('Time (min)')
    ylabel('Runt protein conc. (AU)')
    %legend('male','female','Location','NorthWest')
    StandardFigure(gcf,gca)

    % Save the plots
    FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\BGsubtracted_by_Cyto\MixedSex';
    saveas(gcf,[FigPath,filesep,'Runt_Time_NC13-14_Mixed @ AP =',num2str(APpos),'%','.tif']); 
    saveas(gcf,[FigPath,filesep,'Runt_Time_NC13-14_Mixed @ AP =',num2str(APpos),'%','.pdf']); 
end
%% Runt conc. over Time (male & female) - NC14 only
for APbin = 9:25 % Covers 20%-60% of the embryo
    clf
    APpos = (APbin - 1) * 2.5;

    hold on
    errorbar(Time(NC14:end) - Time(NC14),...
                Averaged_Fluo(NC14:end,APbin),...
                SE_Fluo(NC14:end,APbin))

    % Vertical bar for NC demonstration
    %xline(NC14,'--k')
    % xlim([0 max(Time_female)])
    xlim([0 30])
    ylim([100 max(Averaged_Fluo(NC14:NC14+30,APbin)) + 50])
    title('Runt protein concentration over Time')
    xlabel('Time (min)')
    ylabel('Runt protein conc. (AU)')
    %legend('male','female','Location','NorthWest')
    StandardFigure(gcf,gca)

    % Save the plots
    FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Input dynamics\BGsubtracted_by_Cyto\MixedSex';
    saveas(gcf,[FigPath,filesep,'Runt_Time_NC14_Mixed @ AP =',num2str(APpos),'%','.tif']); 
    saveas(gcf,[FigPath,filesep,'Runt_Time_NC14_Mixed @ AP =',num2str(APpos),'%','.pdf']); 
end

end