function Plot_OpposingGradient_MS2Data
%% Description : Code by Paul and Yang Joon to plot the MS2 spot fluorescence traces at differents AP

%r0Data = LoadMS2Sets('r0')
%r0Old = load('r0.mat');
% r0=load('r0new.mat');
% r1=load('r1new.mat');
% r2=load('r2new.mat');
% r3=load('r3new.mat');

% r0=load('r0_From_NC12.mat');
% r1=load('r1_From_NC12.mat');
% r2=load('r2_From_NC12.mat');
% r3=load('r3_From_NC12.mat');

FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\TxnOutput_sexed';

r0=load([FilePath,filesep,'r0-new-female.mat']);
r1=load([FilePath,filesep,'r1-new-female.mat']);
r2=load([FilePath,filesep,'r2-new-female.mat']);
r3=load([FilePath,filesep,'r3-new-female.mat']);

%Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt-1min-200Hz-Female-Averaged.mat');
%Bcd = load('Bcd-Averaged.mat');
%% Plot individual datasets with Avearaged (over multiple embryos) to check
% This is now the part of the AverageDatasets_FineTuning.m
%% Ploting the time traces of MS2 spots,Bcd, and Runt
%Choose an AP position

AP_Pos = 15;
NC = 13;

if NC==12
    r0Duration = r0.nc12:r0.nc13;
    r1Duration = r1.nc12:r1.nc13;
    r2Duration = r2.nc12:r2.nc13;
    r3Duration = r3.nc12:r3.nc13;
    %RuntDuration = Runt.nc12:Runt.nc13;
    %BcdDuration = Bcd.nc12:Bcd.nc13;
%     r0Delay = 0;
%     r1Delay = 0;
%     r2Delay = 0;
%     r3Delay = 0;
elseif NC==13
    r0Duration = r0.nc13:r0.nc14;
    r1Duration = r1.nc13:r1.nc14;
    r2Duration = r2.nc13:r2.nc14;
    r3Duration = r3.nc13:r3.nc14;
    %RuntDuration = Runt.nc13:Runt.nc14;
    %BcdDuration = Bcd.nc13:Bcd.nc14;
end
    
% Scaling Bicoid and Runt
RuntScale = 15;
BcdScale = 3;
%BcdBG = min(BcdScale*Bcd.MeanVectorAP(BcdDuration,AP_Pos)); % Background fluo
%t=ElapsedTimeTotal;
%p=plot(t,r0MeanVectorAP(:,AP_Pos),t,r1MeanVectorAP(:,AP_Pos),t,r2MeanVectorAP(:,AP_Pos),t,r3MeanVectorAP(:,AP_Pos),'LineWidth',1);
hold on
errorbar(r0.ElapsedTime(r0Duration),r0.MeanVectorAP(r0Duration,AP_Pos),...
            r0.SDVectorAP(r0Duration,AP_Pos)./sqrt(r0.NParticlesAP(r0Duration,AP_Pos)))
errorbar(r1.ElapsedTime(r1Duration)+1,r1.MeanVectorAP(r1Duration,AP_Pos),...
            r1.SDVectorAP(r1Duration,AP_Pos)./sqrt(r1.NParticlesAP(r1Duration,AP_Pos)))
errorbar(r2.ElapsedTime(r2Duration)+1,r2.MeanVectorAP(r2Duration,AP_Pos),...
            r2.SDVectorAP(r2Duration,AP_Pos)./sqrt(r2.NParticlesAP(r2Duration,AP_Pos)))
errorbar(r3.ElapsedTime(r3Duration),r3.MeanVectorAP(r3Duration,AP_Pos),...
            r3.SDVectorAP(r3Duration,AP_Pos)./sqrt(r3.NParticlesAP(r3Duration,AP_Pos)))
% % Runt
% errorbar(Runt.ElapsedTime(RuntDuration)-Runt.ElapsedTime(RuntDuration(1)),...
%             RuntScale*Runt.MeanVectorAP(RuntDuration,AP_Pos) - min(RuntScale*Runt.MeanVectorAP(RuntDuration,AP_Pos)),...
%             RuntScale*Runt.SDVectorAP(RuntDuration,AP_Pos)./sqrt(Runt.NParticlesAP(RuntDuration,AP_Pos)))
% % Bcd
% errorbar(Bcd.ElapsedTime(BcdDuration) - Bcd.ElapsedTime(BcdDuration(1)),...
%             BcdScale*Bcd.MeanVectorAP(BcdDuration,AP_Pos) - BcdBG,...
%             BcdScale*Bcd.SDVectorAP(BcdDuration,AP_Pos)./sqrt(Bcd.NParticlesAP(BcdDuration,AP_Pos)))

title({['Mean fluorescence during nc13'];['at AP = ',num2str((AP_Pos-1)*2.5),'%']})
xlabel('Time (min)')
ylabel('Mean fluorescence (AU)')
xlim([0 20])
%ylim([0 1400])
legend('r0','r1','r2','r3')%,'Runt','Bcd')

StandardFigure(gcf,gca)
%standardizeFigure(gca,legend,'purple','magenta','yellow','red','blue','green')
%saveas(['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Figures-OpposingGradients\Figures_for_Paper','Mean fluorescence',num2str((AP_Pos-1)*2.5),'%.fig'])
%saveas(gcf,'AP20.tif')

%% Plot only the inputs
% Runt
hold on
errorbar(Runt.ElapsedTime(RuntDuration)-Runt.ElapsedTime(RuntDuration(1)),...
            RuntScale*Runt.MeanVectorAP(RuntDuration,AP_Pos) - min(RuntScale*Runt.MeanVectorAP(RuntDuration,AP_Pos)),...
            RuntScale*Runt.SDVectorAP(RuntDuration,AP_Pos)./sqrt(Runt.NParticlesAP(RuntDuration,AP_Pos)))
% Bcd
errorbar(Bcd.ElapsedTime(BcdDuration) - Bcd.ElapsedTime(BcdDuration(1)),...
            BcdScale*Bcd.MeanVectorAP(BcdDuration,AP_Pos) - BcdBG,...
            BcdScale*Bcd.SDVectorAP(BcdDuration,AP_Pos)./sqrt(Bcd.NParticlesAP(BcdDuration,AP_Pos)))
        

title({['Input TF concentration during nc13'];['at AP = ',num2str((AP_Pos-1)*2.5),'%']})
xlabel('Time (min)')
ylabel('Input TF concentration (AU)')
xlim([0 20])
ylim([0 1400])
legend('Runt','Bcd')

standardizeFigure(gca,legend,[])

%% Plot only the output
hold on
errorbar(r0.ElapsedTime(r0Duration),r0.MeanVectorAP(r0Duration,AP_Pos),...
            r0.SDVectorAP(r0Duration,AP_Pos)./sqrt(r0.NParticlesAP(r0Duration,AP_Pos)))
errorbar(r1.ElapsedTime(r1Duration),r1.MeanVectorAP(r1Duration,AP_Pos),...
            r1.SDVectorAP(r1Duration,AP_Pos)./sqrt(r1.NParticlesAP(r1Duration,AP_Pos)))
errorbar(r2.ElapsedTime(r2Duration),r2.MeanVectorAP(r2Duration,AP_Pos),...
            r2.SDVectorAP(r2Duration,AP_Pos)./sqrt(r2.NParticlesAP(r2Duration,AP_Pos)))
errorbar(r3.ElapsedTime(r3Duration),r3.MeanVectorAP(r3Duration,AP_Pos),...
            r3.SDVectorAP(r3Duration,AP_Pos)./sqrt(r3.NParticlesAP(r3Duration,AP_Pos)))
        
title({['Mean spot fluorescence during nc13'];['at AP = ',num2str((AP_Pos-1)*2.5),'%']})
xlabel('Time (min)')
ylabel('Mean spot fluorescence (AU)')
xlim([0 20])
ylim([0 850])
legend('r0','r1','r2','r3')

standardizeFigure(gca,legend,[])
%% Plot the Fraction ON
% Here, I will plot the FractionON that are calculated as the
% Sum of all ON Nuclei / Sum of all Nuclei (sum over multiple embryos)
for NC=12:14
    NCIndex = NC-11;
    figure(NCIndex)
    hold on
    errorbar(0:0.025:1,r0.FractionON_Average(:,NCIndex),r0.FractionON_Average_Error(:,NCIndex))
    errorbar(0:0.025:1,r1.FractionON_Average(:,NCIndex),r1.FractionON_Average_Error(:,NCIndex))
    errorbar(0:0.025:1,r2.FractionON_Average(:,NCIndex),r2.FractionON_Average_Error(:,NCIndex))
    errorbar(0:0.025:1,r3.FractionON_Average(:,NCIndex),r3.FractionON_Average_Error(:,NCIndex))
    xlim([0 0.8])
    title(['Fraction of Active Nuclei ','in NC=',num2str(NC)])
    xlabel('AP')
    ylabel('Fraction of Active Nuclei')
    legend('r0','r1','r2','r3')
    standardizeFigure_YJK(gca,legend,'blue')
end

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
end