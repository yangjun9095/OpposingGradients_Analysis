% generate plots of MS2 from a single embryo
% Last updated : 05/27/2021
% Note that the mean spot fluorescence trace is over competent nuclei, not
% instantaneously on.

%projectName = '2019-07-19-hbP2-r0_MS2V5-lacZ-4';
% projectName = '2019-07-25-hbP2-r0_MS2V5-lacZ-5';
projectName = '2020-08-01-run3-hbP2-r2far-MS2V5-lacZ-MCP-GFP-1';

% load the dataset
DataPath = 'S:\YangJoon\Dropbox\OpposingGradient';
cp = load([DataPath, filesep, projectName, filesep, 'CompiledParticles.mat']);

% extract useful fields
Time = cp.ElapsedTime;
nc13 = cp.nc13;
nc14 = cp.nc14;

Fluo_mean = cp.MeanVectorAP{1,1};
Fluo_SD = cp.SDVectorAP{1,1};
NParticles = cp.NParticlesAP{1,1};
Fluo_SE = Fluo_SD./NParticles;
OnRatio = cell2mat(cp.OnRatioAP);



%% generate plots

% pick an AP bin and time frames
APbin = 11; % 25%
% define the time window (nc14 or nc13)
tWindow = nc14-2:length(Time);
errorbar(Time(tWindow) - Time(nc14)+1.3,...
            Fluo_mean(tWindow,APbin).*OnRatio(tWindow,APbin),...
            Fluo_SE(tWindow,APbin).*OnRatio(tWindow,APbin))
box on
xlabel(' time into nc14 (min)')
ylabel('RNAP number (AU)')

xlim([0 40])
xticks([0 5 10 15 20 25 30 35 40])
ylim([0 1200])

StandardFigure(gcf,gca)

% save the plot
figPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\InitialSlopes\trace_example_101';
saveas(gcf,[figPath,filesep,'101_MS2_TimeTrace_25%','.tif']); 
saveas(gcf,[figPath,filesep,'101_MS2_TimeTrace_25%','.pdf']); 