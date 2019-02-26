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

r0=load('r0-new.mat');
r1=load('r1-new.mat');
r2=load('r2-new.mat');
r3=load('r3-new.mat');

Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt-1min-200Hz-Female-Averaged.mat');
Bcd = load('Bcd-Averaged.mat');
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
    RuntDuration = Runt.nc12:Runt.nc13;
    BcdDuration = Bcd.nc12:Bcd.nc13;
%     r0Delay = 0;
%     r1Delay = 0;
%     r2Delay = 0;
%     r3Delay = 0;
elseif NC==13
    r0Duration = r0.nc13:r0.nc14;
    r1Duration = r1.nc13:r1.nc14;
    r2Duration = r2.nc13:r2.nc14;
    r3Duration = r3.nc13:r3.nc14;
    RuntDuration = Runt.nc13:Runt.nc14;
    BcdDuration = Bcd.nc13:Bcd.nc14;
end
    
% Scaling Bicoid and Runt
RuntScale = 15;
BcdScale = 3;
BcdBG = min(BcdScale*Bcd.MeanVectorAP(BcdDuration,AP_Pos)); % Background fluo
%t=ElapsedTimeTotal;
%p=plot(t,r0MeanVectorAP(:,AP_Pos),t,r1MeanVectorAP(:,AP_Pos),t,r2MeanVectorAP(:,AP_Pos),t,r3MeanVectorAP(:,AP_Pos),'LineWidth',1);
hold on
errorbar(r0.ElapsedTime(r0Duration),r0.MeanVectorAP(r0Duration,AP_Pos),...
            r0.SDVectorAP(r0Duration,AP_Pos)./sqrt(r0.NParticlesAP(r0Duration,AP_Pos)))
errorbar(r1.ElapsedTime(r1Duration),r1.MeanVectorAP(r1Duration,AP_Pos),...
            r1.SDVectorAP(r1Duration,AP_Pos)./sqrt(r1.NParticlesAP(r1Duration,AP_Pos)))
errorbar(r2.ElapsedTime(r2Duration),r2.MeanVectorAP(r2Duration,AP_Pos),...
            r2.SDVectorAP(r2Duration,AP_Pos)./sqrt(r2.NParticlesAP(r2Duration,AP_Pos)))
errorbar(r3.ElapsedTime(r3Duration),r3.MeanVectorAP(r3Duration,AP_Pos),...
            r3.SDVectorAP(r3Duration,AP_Pos)./sqrt(r3.NParticlesAP(r3Duration,AP_Pos)))
% Runt
errorbar(Runt.ElapsedTime(RuntDuration)-Runt.ElapsedTime(RuntDuration(1)),...
            RuntScale*Runt.MeanVectorAP(RuntDuration,AP_Pos) - min(RuntScale*Runt.MeanVectorAP(RuntDuration,AP_Pos)),...
            RuntScale*Runt.SDVectorAP(RuntDuration,AP_Pos)./sqrt(Runt.NParticlesAP(RuntDuration,AP_Pos)))
% Bcd
errorbar(Bcd.ElapsedTime(BcdDuration) - Bcd.ElapsedTime(BcdDuration(1)),...
            BcdScale*Bcd.MeanVectorAP(BcdDuration,AP_Pos) - BcdBG,...
            BcdScale*Bcd.SDVectorAP(BcdDuration,AP_Pos)./sqrt(Bcd.NParticlesAP(BcdDuration,AP_Pos)))

title({['Mean fluorescence during nc13'];['at AP = ',num2str((AP_Pos-1)*2.5),'%']})
xlabel('Time (min)')
ylabel('Mean fluorescence (AU)')
xlim([0 20])
ylim([0 1400])
legend('r0','r1','r2','r3','Runt','Bcd')

standardizeFigure(gca,legend,'purple','magenta','yellow','red','blue','green')
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

end