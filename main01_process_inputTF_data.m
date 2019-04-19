function main01_process_inputTF_data
% This script is for 
% 1) processing the input TF datasets, meaning compiling, averaging, etc.
% 2) plotting the input dynamics of Bicoid and Runt over time and space.

%% Averaging multiple datasets
% Let's use AverageDatasets_NuclearProtein.m script and use DataStatus.xlsx
% INPUT : datasets saved in Dropbox folder, 
% OUTPUT : Check the file name, .mat type

% Averaging Runt datasets (male/female)
%AverageDatasets_NuclearProtein('Runt-1min-200Hz-Male','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
%AverageDatasets_NuclearProtein('Runt-1min-200Hz-Female','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');

% Averaging Bcd and Runt datasets (this is using old Runt datasets, female,
% hets for JB3-Runt/wt Runt, just to see the dynamics)
AverageDatasets_NuclearProtein('Bcd','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');
AverageDatasets_NuclearProtein('Runt','NC',13,'savePath','E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData');

%% Load datasets (Bcd and Runt) - the Runt datasets shoudl be replaced with male/female datasets.
BcdData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Bcd-Averaged.mat')

RuntData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-Averaged.mat')
RuntData_female = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\Runt-1min-200Hz-Female-Averaged.mat')
%% Extract useful fields (Bcd)
BcdTime = BcdData.ElapsedTime;
BcdFluo = BcdData.MeanVectorAP;
BcdFluoSD = BcdData.SDVectorAP;
BcdFluoNNuclei = BcdData.NParticlesAP;
BcdFluoSE = BcdFluoSD./BcdFluoNNuclei;

BcdNC12 = BcdData.nc12;
BcdNC13 = BcdData.nc13;
BcdNC14 = BcdData.nc14;

%% Extract useful fields (Runt)
RuntData = RuntData_female;
% female
RuntTime = RuntData.ElapsedTime;
RuntFluo = RuntData.MeanVectorAP;
RuntFluoSD = RuntData.SDVectorAP;
RuntFluoNNuclei = RuntData.NParticlesAP;
RuntFluoSE = RuntFluoSD./RuntFluoNNuclei;

RuntNC12 = RuntData.nc12;
RuntNC13 = RuntData.nc13;
RuntNC14 = RuntData.nc14;


%% Plot over AP
BcdLength = length(BcdNC13:length(BcdTime));
RuntLength = length(RuntNC13:length(RuntTime));

tLength = min(BcdLength,RuntLength);

% Plot
hold on
for i=[RuntNC12+10,RuntNC12+20,...
        RuntNC13+10,RuntNC13+30,...
        RuntNC14+10,RuntNC14+40,RuntNC14+70,RuntNC14+90];
    errorbar(0:0.025:1,RuntFluo(i,:),RuntFluoSD(i,:))
    title('Runt concentration')
    xlabel('AP')
    ylabel('Runt concentration (AU)')
    ylim([0 500])
    pause
end
legend('NC12 early','NC12 late','NC13 early','NC13 late','NC14 early','NC14 mid','NC14 late')
%% Plot Bcd and Runt over Time
AP = 17;

BcdBG = min(BcdFluo(:,AP)); % This needs to be measured carefully later, for the autofluorescence

hold on
errorbar(BcdTime(BcdNC13:3:tLength)-BcdTime(BcdNC13),BcdFluo(BcdNC13:3:tLength,AP)-BcdBG,BcdFluoSD(BcdNC13:3:tLength,AP))
errorbar(RuntTime(RuntNC13:3:tLength),RuntFluo(RuntNC13:3:tLength,AP),RuntFluoSD(RuntNC13:3:tLength,AP))

title(['Transcription factor concentration',' @ AP = ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Protein concentration (AU)')
legend('Bcd','Runt')
standardizeFigure(gca,legend,[])

end