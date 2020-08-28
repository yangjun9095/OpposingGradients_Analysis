% generate plots of MS2 from a single embryo
% Last updated : 08/20/2020
% 

%projectName = '2019-07-19-hbP2-r0_MS2V5-lacZ-4';
projectName = '2019-07-25-hbP2-r0_MS2V5-lacZ-5';

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



%% generate plots

% pick an AP bin and time frames
APbin = 11; % 25%
% define the time window (nc14 or nc13)
tWindow = nc14:length(Time);
errorbar(Time(tWindow) - Time(nc14),...
            Fluo_mean(tWindow,APbin),...
            Fluo_SE(tWindow,APbin))
box on
xlabel(' time into nc14 (min)')
ylabel('mean spot fluorescence (AU)')

StandardFigure(gcf,gca)

% save the plot
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\OpposingGradientsFigures\Data\Ms2_TimeTraces_example';
saveas(gcf,[figPath,filesep,'000_new_MS2_TimeTrace_25%','.tif']); 
saveas(gcf,[figPath,filesep,'000_new_MS2_TimeTrace_25%','.pdf']); 