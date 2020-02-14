function main02_19_Chek_InitialSlope_individualEmbryos
%% DESCRIPTION
% In this script, I'll double check the initial slopes (fitted) for
% individual embryos by plotting the RateFit with their own SDRateFit
% (error estimated by Confidence Intervals).

%% Predecing steps :
% main02_08_plot_InitialSlopes_AllConstructs.m
% Define the datasets to process
load(['E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData',filesep,'InitialSlopes_ONnuclei_AllConstructs.mat'],...
    'initialSlopes_ONnuclei','-v7.3');

%% Color sets
% For the averaged one : use blue
% For the individual ones : use red-ish colors with some gradation.

red_color{1} = [117,57,73]/255;
red_color{2} = [145,39,30]/255;
red_color{3} = [162,77,55]/255;
red_color{4} = [154,36,123]/255;
red_color{5} = [129,41,140]/255;
red_color{6} = [80,20,43]/255;
red_color{7} = [219,194,197]/255;

%% Plot the initial slope (over ON nuclei) averaged over embryos vs individual embryos
% This is to check how much embryo-to-embryo variability there is.
APaxis = 0:0.025:1;
legendTags = {'embryo1','embryo2','embryo3','embryo4','embryo5','embryo6'}
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\Average_Individual'

for i=2:length(initialSlopes_ONnuclei) % for each data type
    datasetName = initialSlopes_ONnuclei{i,1};
    fittedRate = initialSlopes_ONnuclei{i,2};
    fittedRate = squeeze(fittedRate(:,3,:)); % Extract the NC14
    
    fittedRateSD = initialSlopes_ONnuclei{i,3};
    fittedRateSD = squeeze(fittedRateSD(:,3,:)); % Extract the NC14
    
    average_fittedRate = initialSlopes_ONnuclei{i,5};
    average_fittedRate = squeeze(average_fittedRate(:,3)); % Extract the NC14
    
    SEM_fittedRate = initialSlopes_ONnuclei{i,6};
    SEM_fittedRate = squeeze(SEM_fittedRate(:,3)); % Extract the NC14
    
    % plot averaged profile
    figure(i)
    hold on
    errorbar(APaxis, average_fittedRate, SEM_fittedRate)
    % plot individual embryos
    for j=1:length(fittedRate(1,:)) % indexing embryos
        errorbar(APaxis, fittedRate(:,j), fittedRateSD(:,j),'Color',red_color{j})
    end
    
    xlim([0.15 0.6])
    title(['initial slope (ON nuclei)-',datasetName])
    xlabel('AP axis (EL)')
    ylabel('initial slope (AU/min)')
    legend('averaged')
    StandardFigure(figure(i), figure(i).CurrentAxes)
    saveas(figure(i),[FigPath,filesep, 'InitialRates_Averaged_Individual',...
                        '_NC14_ONnuclei_',datasetName , '.pdf']); 
    saveas(figure(i),[FigPath,filesep, 'InitialRates_Averaged_Individual',...
                        '_NC14_ONnuclei_',datasetName , '.tif']); 
end