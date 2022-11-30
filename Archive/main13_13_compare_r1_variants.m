%% Script to compare r1 expression at different positions (r1 vs r1-close)
function main13_13_compare_r1_variants
%% Process the datasets
% Use AverageDatasets.m, and AccumulatedmRNA.m for averaging, and
% accumulating the mRNA.
AverageDatasets('r1-close','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData')

% Accumulate mRNA
AccumulatedmRNA('r1-close',2)
%% Load the datasets
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
Data_r1 = load([filePath,filesep,'r1-new-female.mat']); % female
Data_r1_close = load([filePath,filesep,'r1-close.mat']); % mixed, but mostly female.

% AccumulatedmRNA_r1 = load([filePath,filesep,'AccumulatedmRNA_r1-new-female.mat']);
% AccumulatedmRNA_r1_close = load([filePath,filesep,'AccumulatedmRNA_r1-close.mat'])
%% Extract useful field
% r1
Time_r1 = Data_r1.ElapsedTime;
NC13_r1 = Data_r1.nc13;
NC14_r1 = Data_r1.nc14;
MeanVectorAP_r1 = Data_r1.MeanVectorAP;
SDVectorAP_r1 = Data_r1.SDVectorAP;
SEVectorAP_r1 = Data_r1.SEVectorAP;
NParticlesAP_r1 = Data_r1.NParticlesAP;

% r1_close
Time_r1_close = Data_r1_close.ElapsedTime;
NC13_r1_close = Data_r1_close.nc13;
NC14_r1_close = Data_r1_close.nc14;
MeanVectorAP_r1_close = Data_r1_close.MeanVectorAP;
SDVectorAP_r1_close = Data_r1_close.SDVectorAP;
SEVectorAP_r1_close = Data_r1_close.SEVectorAP;
NParticlesAP_r1_close = Data_r1_close.NParticlesAP;

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.purple; colorDict.green; colorDict.brown; colorDict.red; colorDict.brown]; % 4 embryos max. it could be extended easily
%% 1) Plot for mean spot fluorescence _ individual embryos
APbin = 1 ; % APbin position

Time = Data_r1_close.ElapsedTime;
MeanVectorAP_individual = Data_r1_close.MeanVectorAP_individual;
SDVectorAP_individual = Data_r1_close.SDVectorAP_individual;

hold on
for i=1:length(MeanVectorAP_individual(1,1,:))
    errorbar(Time, MeanVectorAP_individual(:,APbin,i), SDVectorAP_individual(:,APbin,i))
    pause
end


%% 2) Plot for mean spot fluorescence
APbin = 15; % APbin position

% NC13
Range_r1 = NC13_r1:NC14_r1;
Range_r1_close = NC13_r1_close:NC14_r1_close;
tDelay = zeros(1,2);

% % NC14
% Range_r1 = NC14_r1:length(Time_r1);
% Range_r1_close = NC14_r1_close:length(Time_r1_close);
% tDelay = [Time_r1(NC14_r1), Time_r1_close(NC14_r1_close)];

MS2TraceFig = figure;
hold on    
errorbar(Time_r1(Range_r1)-tDelay(1), MeanVectorAP_r1(Range_r1,APbin),...
            SEVectorAP_r1(Range_r1,APbin),'Color',ColorChoice(1,:))
        
errorbar(Time_r1_close(Range_r1_close)-tDelay(2), MeanVectorAP_r1_close(Range_r1_close,APbin),...
            SEVectorAP_r1_close(Range_r1_close,APbin),'Color',ColorChoice(2,:))

% xlim, ylim
%xlim([0 20])
ylim([0 max(MeanVectorAP_r1_close(Range_r1_close,APbin)) + 200])
title(['Mean spot fluorescence over time @ AP = ',num2str((APbin-1)*2.5),'%'])
xlabel('Time into NC13 (min)')
ylabel('Mean spot fluorescence (AU)')
legend('r1','r1-close')
StandardFigure(MS2TraceFig,MS2TraceFig.CurrentAxes)

% Save Figure
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\hbP2-r0123-Averaged_MS2_traces_multipleEmbryos';
% saveas(MS2TraceFig,[FigPath,filesep, 'Mean_MS2_traces at AP=', num2str((APbin-1)*2.5),'%' , '_NC14' , '.tif']); 
% saveas(MS2TraceFig,[FigPath,filesep, 'Mean_MS2_traces at AP=', num2str((APbin-1)*2.5),'%' , '_NC14' , '.pdf']); 

%% 2. Accumulated mRNA
%% Alternative calculation for the Accumulated mRNA
% % Since my script, AverageDatasets didn't take into account of the
% % APbinArea for calculating the total mRNA
% % I'll try IntegratemRNA.m script for this.
% 
% 1) Load the datasets
r1Data = LoadMS2Sets('r1-new-female','dontCompare')
r1_close_Data = LoadMS2Sets('r1-close','dontCompare')

% 2) IntegratemRNA

[TotalProd_r1,TotalProdError_r1,TotalProdN_r1,...
    MeanTotalProd_r1,SDTotalProd_r1,SETotalProd_r1]=IntegratemRNA(r1Data,1,2)

[TotalProd_r1_close,TotalProdError_r1_close,TotalProdN_r1_close,...
    MeanTotalProd_r1_close,SDTotalProd_r1_close,SETotalProd_r1_close]=IntegratemRNA(r1_close_Data,1,2)
%% Plot the Accumulated mRNA
APaxis = 0:0.025:1;
% Figure Path
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Transcription-Output\hbP2-r0123-AccumulatedmRNA-females';

% NC13
NC= 13;
AccumulatedmRNA_NC13_figure = figure
hold on
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r1_close(:,NC), SETotalProd_r1_close(:,NC))
xlim([0.2 0.5])

title('Accumulated mRNA over AP @ NC13')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r1','r1-close')

StandardFigure(AccumulatedmRNA_NC13_figure,AccumulatedmRNA_NC13_figure.CurrentAxes)
%saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.tif'])
%saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.pdf'])

% NC14
NC = 14;
AccumulatedmRNA_NC14_figure = figure
hold on
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r1_close(:,NC), SETotalProd_r1_close(:,NC))
xlim([0.2 0.5])

title('Accumulated mRNA over AP @ NC14')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r1','r1-close')

StandardFigure(AccumulatedmRNA_NC14_figure,AccumulatedmRNA_NC14_figure.CurrentAxes)
%saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.tif'])
%saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.pdf'])
end