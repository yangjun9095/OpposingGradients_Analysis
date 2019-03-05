function main01_process_inputTF_data
% This script is for plotting the input dynamics of Bicoid and Runt over
% time and space.

%% Averaging multiple datasets
% Let's use AverageDatasets_NuclearProtein.m script and use DataStatus.xlsx
% INPUT : datasets saved in Dropbox folder, 
% OUTPUT : 


%% Load datasets (Bcd and Runt)
BcdData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\Bcd-Averaged.mat')
RuntData = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\2018-08-19-RuntN-JB3-6-vasa-eGFP1\CompiledNuclei.mat')

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
    title('Runt concentration (BG subtracted with NoNB dataset)')
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