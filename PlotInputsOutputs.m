function PlotInputsOutputs
%% Plotting two inputs (Bcd and Runt), and MS2 output over time
% First, I will pick several AP bins, then plot altogether.

%% Load the datasets
Bcd = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\eGFP-Bcd-From-Liz-Jonathan\BcdGFPAnt.mat')
Bcd = Bcd.DataBcd;
Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-04-29-Runt-JB3-MCP-mCherry-vasa-eGFP1-edited-2.5%APbins\CompiledNuclei.mat')

r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0.mat');
r1Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1.mat');
r2Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2.mat');
r3Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3.mat');
%% Truncate the time series, and synchronize to the beginning of nc13
BcdTime = Bcd.ElapsedTime(Bcd.nc13:end)-Bcd.ElapsedTime(Bcd.nc13);
BcdFluo = Bcd.MeanVectorAP;
BcdFluo = BcdFluo(Bcd.nc13:end,:);
BcdFluo(isnan(BcdFluo)) = 0;
BcdFluoSE = Bcd.SDVectorAP./sqrt(Bcd.NParticlesAP);
BcdFluoSE = BcdFluoSE(Bcd.nc13:end,:)


RuntTime = Runt.ElapsedTime(Runt.nc13:end)-Runt.ElapsedTime(Runt.nc13);
RuntFluo = Runt.MeanVectorAP;
RuntFluo = RuntFluo(Runt.nc13:end,:);
RuntFluo = RuntFluo - 95; % Subtracting the Background for now.
RuntFluo(isnan(RuntFluo)) = 0;
RuntFluoSE = Runt.SDVectorAP./sqrt(Runt.NParticlesAP);
RuntFluoSE = RuntFluoSE(Runt.nc13:end,:);

%% r0
r0Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r0.mat')
Time_r0 = r0Data.ElapsedTime;
MeanVectorAP_r0 = r0Data.MeanVectorAP;
SDVectorAP_r0 = r0Data.SDVectorAP;
SEVectorAP_r0 = r0Data.SEVectorAP;
NParticlesAP_r0 = r0Data.NParticlesAP;
AccumulatedmRNA_r0 = r0Data.AccumulatedmRNA;
SDAccumulatedmRNA_r0 = r0Data.AccumulatedmRNA_SD;

%% r1
r1Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1.mat')
Time_r1 = r1Data.ElapsedTime;
MeanVectorAP_r1 = r1Data.MeanVectorAP;
SDVectorAP_r1 = r1Data.SDVectorAP;
SEVectorAP_r1 = r1Data.SEVectorAP;
NParticlesAP_r1 = r1Data.NParticlesAP;
AccumulatedmRNA_r1 = r1Data.AccumulatedmRNA;
SDAccumulatedmRNA_r1 = r1Data.AccumulatedmRNA_SD;

%% r2
r2Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r2.mat')
Time_r2 = r2Data.ElapsedTime;
MeanVectorAP_r2 = r2Data.MeanVectorAP;
SDVectorAP_r2 = r2Data.SDVectorAP;
SEVectorAP_r2 = r2Data.SEVectorAP;
NParticlesAP_r2 = r2Data.NParticlesAP;
AccumulatedmRNA_r2 = r2Data.AccumulatedmRNA;
SDAccumulatedmRNA_r2 = r2Data.AccumulatedmRNA_SD;

%% r3
r3Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r3.mat')
Time_r3 = r3Data.ElapsedTime;
MeanVectorAP_r3 = r3Data.MeanVectorAP;
SDVectorAP_r3 = r3Data.SDVectorAP;
SEVectorAP_r3 = r3Data.SEVectorAP;
NParticlesAP_r3 = r3Data.NParticlesAP;
AccumulatedmRNA_r3 = r3Data.AccumulatedmRNA;
SDAccumulatedmRNA_r3 = r3Data.AccumulatedmRNA_SD;
%% Plot
AP = 9; %9(20%), 17(40%), 25(60%), 21(50%)
BcdScale = 5;
RuntScale = 5;
hold on
errorbar(BcdTime,BcdScale*BcdFluo(:,AP),BcdScale*BcdFluoSE(:,AP))
errorbar(RuntTime,RuntScale*RuntFluo(:,AP),RuntScale*RuntFluoSE(:,AP))

errorbar(Time_r0,MeanVectorAP_r0(:,AP),SEVectorAP_r0(:,AP))
%errorbar(Time_r1,MeanVectorAP_r1(:,AP),SEVectorAP_r1(:,AP))
%errorbar(Time_r2,MeanVectorAP_r2(:,AP),SEVectorAP_r2(:,AP))
errorbar(Time_r3,MeanVectorAP_r3(:,AP),SEVectorAP_r3(:,AP))

title(['Protein inputs and MS2 output over time',' @ AP = ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Fluorescence (AU)')
legend('Bicoid','Runt','r0','r1','r2','r3')

% standardize figure
standardizeFigure(gca,legend,'b','r','yellow','cyan','magenta','lightblue')
end