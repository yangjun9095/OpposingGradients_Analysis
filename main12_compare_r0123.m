function main12_compare_r0123(varargin)
%% DESCRIPTION
% This script is to compare 
% 1) Mean spot fluo
% 2) integrated mRNA profiles (over AP, at different NC ranges, etc.) 
% between males and females of the same synthetic enhancers.

% The plan is to use this script to compare r0, r1, r2, and r3 

%% Load datasets
% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r3-new-female

FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\TxnOutput_sexed';

% This is because our new datasets(r0-new-male, r0-new-female lack NC14)
r0Data = load([FilePath, filesep, 'r0']);
%r0Data = load([FilePath, filesep, 'r0']);
%r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\r0.mat')
r1Data = load([FilePath, filesep, 'r1-new-female']);
r2Data = load([FilePath, filesep, 'r2-new-female']);
r3Data = load([FilePath, filesep, 'r3-new-female']);

%% Assign fields to variables

% r0
Time_r0 = r0Data.ElapsedTime;
NC13_r0 = r0Data.nc13;
NC14_r0 = r0Data.nc14;
MeanVectorAP_r0 = r0Data.MeanVectorAP;
SDVectorAP_r0 = r0Data.SDVectorAP;
SEVectorAP_r0 = r0Data.SEVectorAP;
NParticlesAP_r0 = r0Data.NParticlesAP;
AccumulatedmRNA_All_r0 = r0Data.AccumulatedmRNA_FractionON;
AccumulatedmRNA_All_SD_r0 = r0Data.AccumulatedmRNA_FractionON_SD;

% r1
Time_r1 = r1Data.ElapsedTime;
NC13_r1 = r1Data.nc13;
NC14_r1 = r1Data.nc14;
MeanVectorAP_r1 = r1Data.MeanVectorAP;
SDVectorAP_r1 = r1Data.SDVectorAP;
SEVectorAP_r1 = r1Data.SEVectorAP;
NParticlesAP_r1 = r1Data.NParticlesAP;
AccumulatedmRNA_All_r1 = r1Data.AccumulatedmRNA_FractionON;
AccumulatedmRNA_All_SD_r1 = r1Data.AccumulatedmRNA_FractionON_SD;

% r2
Time_r2 = r2Data.ElapsedTime;
NC13_r2 = r2Data.nc13;
NC14_r2 = r2Data.nc14;
MeanVectorAP_r2 = r2Data.MeanVectorAP;
SDVectorAP_r2 = r2Data.SDVectorAP;
SEVectorAP_r2 = r2Data.SEVectorAP;
NParticlesAP_r2 = r2Data.NParticlesAP;
AccumulatedmRNA_All_r2 = r2Data.AccumulatedmRNA_FractionON;
AccumulatedmRNA_All_SD_r2 = r2Data.AccumulatedmRNA_FractionON_SD;

% r3
Time_r3 = r3Data.ElapsedTime;
NC13_r3 = r3Data.nc13;
NC14_r3 = r3Data.nc14;
MeanVectorAP_r3 = r3Data.MeanVectorAP;
SDVectorAP_r3 = r3Data.SDVectorAP;
SEVectorAP_r3 = r3Data.SEVectorAP;
NParticlesAP_r3 = r3Data.NParticlesAP;
AccumulatedmRNA_All_r3 = r3Data.AccumulatedmRNA_FractionON;
AccumulatedmRNA_All_SD_r3 = r3Data.AccumulatedmRNA_FractionON_SD;

%% Plot the MeanVectorAP
APbin = 10; % APbin position

% NC13
Range_r0 = NC13_r0:NC14_r0;
Range_r1 = NC13_r1:NC14_r1;
Range_r2 = NC13_r2:NC14_r2;
Range_r3 = NC13_r3:NC14_r3;
tDelay = zeros(1,4);

% NC14
% Range_r0 = NC14_r0:length(Time_r0);
% Range_r1 = NC14_r1:length(Time_r1);
% Range_r2 = NC14_r2:length(Time_r2);
% Range_r3 = NC14_r3:length(Time_r3);
%tDelay = [Time_r0(NC14_r0), Time_r1(NC14_r1), Time_r2(NC14_r2), Time_r3(NC14_r3)];

hold on
errorbar(Time_r0(Range_r0)-tDelay(1), MeanVectorAP_r0(Range_r0,APbin),...
            SEVectorAP_r0(Range_r0,APbin))
        
errorbar(Time_r1(Range_r1)-tDelay(2), MeanVectorAP_r1(Range_r1,APbin),...
            SEVectorAP_r1(Range_r1,APbin))
        
errorbar(Time_r2(Range_r2)-tDelay(3), MeanVectorAP_r2(Range_r2,APbin),...
            SEVectorAP_r2(Range_r2,APbin))
        
errorbar(Time_r3(Range_r3)-tDelay(4), MeanVectorAP_r3(Range_r3,APbin),...
            SEVectorAP_r3(Range_r3,APbin))
% xlim, ylim
%xlim([0 20])
ylim([0 1000])
title(['Mean spot fluorescence over time @ AP = ',num2str((APbin-1)*2.5),'%'])
xlabel('Time into NC14 (min)')
ylabel('Mean spot fluorescence (AU)')
legend('r0','r1','r2','r3')
StandardFigure(gcf,gca)
% standardizeFigure_YJK(gca,legend,[])

%% Plot the Accumulated mRNA
APaxis = 0:0.025:1;


AccumulatedmRNA_All_r0(isnan(AccumulatedmRNA_All_SD_r0)) = nan;
AccumulatedmRNA_All_r1(isnan(AccumulatedmRNA_All_SD_r1)) = nan;
AccumulatedmRNA_All_r2(isnan(AccumulatedmRNA_All_SD_r2)) = nan;
AccumulatedmRNA_All_r3(isnan(AccumulatedmRNA_All_SD_r3)) = nan;

% NC13
AccumulatedmRNA_NC13_figure = figure
hold on
errorbar(APaxis, AccumulatedmRNA_All_r0(NC14_r0,:), AccumulatedmRNA_All_SD_r0(NC14_r0,:))
errorbar(APaxis, AccumulatedmRNA_All_r1(NC14_r1,:), AccumulatedmRNA_All_SD_r1(NC14_r1,:))
errorbar(APaxis, AccumulatedmRNA_All_r2(NC14_r2,:), AccumulatedmRNA_All_SD_r2(NC14_r2,:))
errorbar(APaxis, AccumulatedmRNA_All_r3(NC14_r3,:), AccumulatedmRNA_All_SD_r3(NC14_r3,:))

title('Accumulated mRNA over AP @ NC13')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3')

AccumulatedmRNA_NC14_figure = figure
hold on
errorbar(APaxis, AccumulatedmRNA_All_r0(end,:), AccumulatedmRNA_All_SD_r0(end,:))
errorbar(APaxis, AccumulatedmRNA_All_r1(end,:), AccumulatedmRNA_All_SD_r1(end,:))
errorbar(APaxis, AccumulatedmRNA_All_r2(end,:), AccumulatedmRNA_All_SD_r2(end,:))
errorbar(APaxis, AccumulatedmRNA_All_r3(end,:), AccumulatedmRNA_All_SD_r3(end,:))

title('Accumulated mRNA over AP @ NC14')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3')

%% Alternative calculation for the Accumulated mRNA
% Since my script, AverageDatasets didn't take into account of the
% APbinArea for calculating the total mRNA
% I'll try IntegratemRNA.m script for this.

% 1) Load the datasets
r0Data = LoadMS2Sets('r0')
r1Data = LoadMS2Sets('r1-new-female')
r2Data = LoadMS2Sets('r2-new-female')
r3Data = LoadMS2Sets('r3-new-female')

% 2) IntegratemRNA
[TotalProd_r0,TotalProdError_r0,TotalProdN_r0,...
    MeanTotalProd_r0,SDTotalProd_r0,SETotalProd_r0]=IntegratemRNA(r0Data,1,2)

[TotalProd_r1,TotalProdError_r1,TotalProdN_r1,...
    MeanTotalProd_r1,SDTotalProd_r1,SETotalProd_r1]=IntegratemRNA(r1Data,1,2)

[TotalProd_r2,TotalProdError_r2,TotalProdN_r2,...
    MeanTotalProd_r2,SDTotalProd_r2,SETotalProd_r2]=IntegratemRNA(r2Data,1,2)

[TotalProd_r3,TotalProdError_r3,TotalProdN_r3,...
    MeanTotalProd_r3,SDTotalProd_r3,SETotalProd_r3]=IntegratemRNA(r3Data,1,2)

%% Plot the Accumulated mRNA
APaxis = 0:0.025:1;
% Figure Path
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Transcription-Output\hbP2-r0123-AccumulatedmRNA-females';
% 
% AccumulatedmRNA_All_r0(isnan(AccumulatedmRNA_All_SD_r0)) = nan;
% AccumulatedmRNA_All_r1(isnan(AccumulatedmRNA_All_SD_r1)) = nan;
% AccumulatedmRNA_All_r2(isnan(AccumulatedmRNA_All_SD_r2)) = nan;
% AccumulatedmRNA_All_r3(isnan(AccumulatedmRNA_All_SD_r3)) = nan;

% NC13
NC= 13;
AccumulatedmRNA_NC13_figure = figure
hold on
errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))
xlim([0.2 0.5])

title('Accumulated mRNA over AP @ NC13')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3')

StandardFigure(AccumulatedmRNA_NC13_figure,AccumulatedmRNA_NC13_figure.CurrentAxes)
%saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.tif'])
%saveas(AccumulatedmRNA_NC13_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC13','_SE','.pdf'])

% NC14
NC = 14;
AccumulatedmRNA_NC14_figure = figure
hold on
errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))

xlim([0.2 0.5])

title('Accumulated mRNA over AP @ NC14')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3')

StandardFigure(AccumulatedmRNA_NC14_figure,AccumulatedmRNA_NC14_figure.CurrentAxes)
%saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.tif'])
%saveas(AccumulatedmRNA_NC14_figure,[FigPath,filesep,'AccumulatedmRNA_r0123',DataType(1:end-1),'_NC14','_SE','.pdf'])

%% Side thing
% Compare the total mRNA profile for diffrent sex
% 1) Load the datasets

r1Data_male = LoadMS2Sets('r1-new-male')
r2Data_male = LoadMS2Sets('r2-new-male')
r3Data_male = LoadMS2Sets('r3-new-male')

% 2) IntegratemRNA

[TotalProd_r1_male,TotalProdError_r1_male,TotalProdN_r1_male,...
    MeanTotalProd_r1_male,SDTotalProd_r1_male,SETotalProd_r1_male]=IntegratemRNA(r1Data_male,1,2)

[TotalProd_r2_male,TotalProdError_r2_male,TotalProdN_r2_male,...
    MeanTotalProd_r2_male,SDTotalProd_r2_male,SETotalProd_r2_male]=IntegratemRNA(r2Data_male,1,2)

[TotalProd_r3_male,TotalProdError_r3_male,TotalProdN_r3_male,...
    MeanTotalProd_r3_male,SDTotalProd_r3_male,SETotalProd_r3_male]=IntegratemRNA(r3Data_male,1,2)

%% Plot for different sexes
NC= 13;
AccumulatedmRNA_NC13_figure = figure
hold on
% females (r0 is mixed)
errorbar(APaxis, MeanTotalProd_r0(:,NC), SETotalProd_r0(:,NC))
errorbar(APaxis, MeanTotalProd_r1(:,NC), SETotalProd_r1(:,NC))
errorbar(APaxis, MeanTotalProd_r2(:,NC), SETotalProd_r2(:,NC))
errorbar(APaxis, MeanTotalProd_r3(:,NC), SETotalProd_r3(:,NC))
% males
errorbar(APaxis, MeanTotalProd_r1_male(:,NC), SETotalProd_r1_male(:,NC))
errorbar(APaxis, MeanTotalProd_r2_male(:,NC), SETotalProd_r2_male(:,NC))
errorbar(APaxis, MeanTotalProd_r3_male(:,NC), SETotalProd_r3_male(:,NC))
xlim([0.2 0.5])

title('Accumulated mRNA over AP @ NC13')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3','r1-male','r2-male','r3-male')
%% Section 2. Mean Spot fluo
%% Mean spot fluorescence between constructs

end