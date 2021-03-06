function main12_compare_r0123(varargin)
%% DESCRIPTION
% This script is to compare the integrated mRNA profiles (over AP, at
% different NC ranges, etc.) between males and females of the same
% synthetic enhancers.

% The plan is to use this script to compare 
% 1) r0, r1, r2, and r3 

%% Load datasets
% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r3-new-female

FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\TxnOutput_sexed';

%r0Data = load([FilePath, filesep, 'r0-new-female']);
r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\r0.mat')
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


end