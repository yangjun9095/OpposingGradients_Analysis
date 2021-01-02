%Plots corner plot of example single-cell MCMC inference

%% Load data
% clear variables
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_inference_ThermoModelV2\MCMCresults';
MCMC = load([FilePath, filesep,'MCMCchain_001.mat']);
MCMCchain = MCMC.MCMCchain;

% extract chain, results, etc. from the saved MCMC result mat file.
% chain = 
% results = 

% extract the n_simu, n_burns
n_simu = results.nsimu;
n_burn = 0.5 * n_simu;

%% Color module
% This is defining the line color
% We have 8 distinct datasets, with or without Runt protein.
% I think selecting 8 distinguishable color sets, then changing the
% brightness by either adding/subtracting white would be a better idea than
% selecting 16 different color sets.

colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.purple = [171,133,172]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

%colorDict.magenta = [208,109,171]/255;
%colorDict.lightBlue = [115,142,193]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.pink = [232,177,157]/255;
colorDict.thickpink = [132,27,69]/255;

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.darkgreen; colorDict.thickpink]; 

% For now, I'll add white (color+[1 1 1])/2 to make thinner color (for the
% Runt nulls)

%% Step0. Set up the directories
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo';
mkdir(FigPath)

%% Raw chain plots

% names
names = results.names;

% number of rows and columns in the subplot
n_rows = 3;
n_cols = 2;

raw = figure;
hold on
box on

for i = 1:length(names)-1
    subplot(n_rows, n_cols, i)
    PlotHandle = plot(chain(n_burn+1:end, i));
    ylabel([names{i}, ' (AU)']);
    xlim([0 100000-n_burn]);
    % set(gca,'XTick',[0 100000],'XScale','log');
    StandardFigure(PlotHandle,gca);
end

% subplot(n_rows,3,1)
% PlotHandle = plot(chain(n_burn+1:end,1));
% ylabel('K_{b} (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);
% 
% subplot(3,3,2)
% PlotHandle = plot(chain(n_burn+1:end,2));
% ylabel('K_{r} (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);
% 
% subplot(3,3,3)
% PlotHandle = plot(chain(n_burn+1:end,3));
% ylabel('\omega_{b} (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);
% 
% subplot(3,3,4)
% PlotHandle = plot(chain(n_burn+1:end,4));
% ylabel('\omega_{bp} (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);
% 
% subplot(3,3,5)
% PlotHandle = plot(chain(n_burn+1:end,5));
% ylabel('\omega_{rp} (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);
% 
% subplot(3,3,6)
% PlotHandle = plot(chain(n_burn+1:end,6));
% ylabel('p (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);
% 
% subplot(3,3,7)
% PlotHandle = plot(chain(n_burn+1:end,7));
% ylabel('R_{max} (AU)');
% xlim([0 100000-n_burn]);
% % set(gca,'XTick',[0 100000],'XScale','log');
% StandardFigure(PlotHandle,gca);

construct = '001'
%Save the plot
saveas(gcf,[FigPath,filesep,'MCMC_raw_chains_', construct ,'.tif']); 
saveas(gcf,[FigPath,filesep,'MCMC_raw_chains_', construct ,'.pdf']);  
%% Raw chain panels
% rawchains = figure;
% mcmcplot(chain, [], results.names, 'chainpanel')

%% Autocorrelation plot
% m = [chain(i).v_chain,chain(i).dwell_chain,chain(i).R_chain(:,1)];
m = [chain(n_burn+1:end,1), chain(n_burn+1:end,2), chain(n_burn+1:end,3), chain(n_burn+1:end,4), chain(n_burn+1:end,5), chain(n_burn+1:end,6)];%, chain(n_burn+1:end,7)];
[C,lags,ESS] = eacorr(m);

acorr = figure;
box on
hold on

PlotHandle = plot(lags,C(:,1),'DisplayName','K_{b}');
PlotHandle(end+1) = plot(lags,C(:,2),'DisplayName','K_{r}');
PlotHandle(end+1) = plot(lags,C(:,3),'DisplayName','\omega_{b}');
PlotHandle(end+1) = plot(lags,C(:,4),'DisplayName','\omega_{bp}');
PlotHandle(end+1) = plot(lags,C(:,5),'DisplayName','\omega_{rp}');
PlotHandle(end+1) = plot(lags,C(:,6),'DisplayName','p');
% PlotHandle(end+1) = plot(lags,C(:,7),'DisplayName','R_{max}');

xlabel('lag');
ylabel('autocorrelation');
ylim([-0.2 1])
set(gca,'XTick',0:50000:n_simu-n_burn);
legend;
StandardFigure(PlotHandle,gca);

saveas(gcf,[FigPath,filesep,'Auto_corr_chains_', construct ,'.tif']); 
saveas(gcf,[FigPath,filesep,'Auto_corr_chains_', construct ,'.pdf']);  
%% Corner plot
n_burn = 0.5*n_steps;
m = [chain(n_burn:end,1), chain(n_burn:end,2), chain(n_burn:end,3), chain(n_burn:end,4), chain(n_burn:end,5), chain(n_burn:end,6)];%, chain(n_burn:end,7)];
corner = figure;
% names = {'Kb','Kr','w_b','w_{bp}','w_{rp}','p','R_{max}'};
ecornerplot(m,'names',names);

saveas(gcf,[FigPath,filesep,'Corner_plot_', construct ,'.tif']); 
saveas(gcf,[FigPath,filesep,'Corner_plot_', construct ,'.pdf']);  
%% Data and inference confidence intervals


%% 
%% Calculate the output range using the mcmcpred
% MCMCdata.xdata
out = mcmcpred(results,chain,[],MCMCdata.xdata, mdl);

%% MCMC prediction plot

% plothandle = mcmcpredplot(out, MCMCdata.ydata, [])

% mcmcpredplot(out);
nn = (size(out.predlims{1}{1},1) + 1) / 2;
plimi = out.predlims{1};
yl = plimi{1}(1,:);
yu = plimi{1}(2*nn-1,:);

% km = mean(MCMCchain);
% ks = std(MCMCchain);
    

yf = plimi{1}(nn,:);
    
yy = yf;
yyl = yl;
yyu = yu;

index1 = 1:11;
index2 = 12:22;

figure
hold on
% Runt null
errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
% shadedErrorBar(APaxis(APbins_fit), yf(1,index1), [yu(1,index1);  yl(1,index1)],'lineProps',{'markerfacecolor',ColorChoice(4,:)})
plot(APaxis(APbins_fit), yf(1,index1), 'Color',ColorChoice(4,:))
plot(APaxis(APbins_fit), yl(1,index1),'Color',(ColorChoice(4,:)+[1 1 1])/2)
plot(APaxis(APbins_fit), yu(1,index1),'Color',(ColorChoice(4,:)+[1 1 1])/2)

% Runt WT
errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
% shadedErrorBar(APaxis(APbins_fit), yf(1,index2), [yu(1,index2);  yl(1,index2)],'lineProps',{'markerfacecolor',ColorChoice(1,:)})
plot(APaxis(APbins_fit), yf(1,index2),'Color',ColorChoice(1,:))
plot(APaxis(APbins_fit), yl(1,index2),'Color',(ColorChoice(1,:)+[1 1 1])/2)
plot(APaxis(APbins_fit), yu(1,index2),'Color',(ColorChoice(1,:)+[1 1 1])/2)                                 

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})

StandardFigure(gcf,gca)

%Save the plot
% saveas(gcf,[FigPath,filesep,'MCMC_fit_direct_001','.tif']); 
% saveas(gcf,[FigPath,filesep,'MCMC_fit_direct_001','.pdf']); 
%% Export graphics
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_inference_ThermoModelV2';

exportgraphics(raw,[FigPath, filesep, 'RawChains.pdf'],'ContentType','vector','BackgroundColor','none');
% exportgraphics(acorr,[FigPath, filesep,'AutoCorrelation.pdf'],'ContentType','vector','BackgroundColor','none');
exportgraphics(corner,[FigPath, filesep,'CornerPlot.pdf'],'ContentType','vector','BackgroundColor','none');