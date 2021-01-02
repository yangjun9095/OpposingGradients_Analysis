%% generate_plots_degeneracy_general_ThermoModel_V2
function gen_plots_degeneracy_general_ThermoModel_V2

%% Descroption
% This script is for generating plots to show the degeneracy in the model.
% Specifically, we want to show how degenerate our general thermodynamic
% model (general_Thermo_model_V2) is, by showing the model fits, and
% parameters for 5-10 cases.

%% Load the fit and params
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo\DirectRepression';

Fits = load([FilePath, filesep, '001_degeneracy_Check.mat'],'FitResult');

Fits_degenerate = Fits.FitResult;

%% pick a couple of fits and params for demonstration of degeneracy
% index = [1, 4, 5, 12, 13, 14, 15, 21, 51];
k=1; % counter
for i = 1:1000
    params_deg(k,:) = Fits_degenerate(i).params_fit;
    fit(k,:) = Fits_degenerate.Ypred;
    fit_SE(k,:) = Fits_degenerate.delta;
    k=k+1;
end

%% generate plots of fits and params
% first, fits

construct = 5; %[001]
% initialize the variables
Rate = [];
Rate_SEM = [];
Rate_null = [];
Rate_null_SEM = [];

% Pull the relevant data from the compiledData
Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% pull the fitted result and error bar (Ypred and delta) from the
% lsqcurvefit fitting.
APrange = 9:19;

hold on
% Runt nulls
errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})
box on
StandardFigure(gcf,gca)

% plot different fits from different sets of parameters
for i=1:1000%length(index)
    plot(APaxis(APrange), fit(i,1:length(APrange)))
    plot(APaxis(APrange), fit(i,1+length(APrange):end))
    pause
end
% shadedErrorBar(APrange, Ypred(1:length(APrange)), delta(1:length(APrange)),'lineProps',{'markerfacecolor',ColorChoice(4,:)})
% shadedErrorBar(APrange, Ypred(length(APrange)+1:end), delta(length(APrange)+1:end),'lineProps',{'markerfacecolor',ColorChoice(1,:)})


%   %Save the plot
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.pdf']); 
end