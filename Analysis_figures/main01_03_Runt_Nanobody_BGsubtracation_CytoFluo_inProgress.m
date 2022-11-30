function main01_03_Runt_Nanobody_BGsubtracation_CytoFluo_inProgress
% This script is for analyzing the LlamaTag(JB3)-NB, eGFP data analysis
% From our fluorescence measurements to get the amount of NB

% The main steps that we need are
% 1) Background fluo subtraction (Note that autofluorescence should be
% considered). This could be done in two different ways. Maybe I should do
% both ways, then see that both show similar results.
% 2) Repeat this for all datasets, both male and females.

%% Load the datasets
RuntMale = LoadMS2Sets('Runt-1min-200Hz-Male');
%NBData = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-05-22-NB-JB3-MCP-mCherry-vasa-eGFP1\CompiledNuclei.mat');

NBData = RuntMale(2);

% Extract some useful fields
NBTime = NBData.ElapsedTime;
NBFluo = NBData.MeanVectorAP;
SDNBFluo = NBData.SDVectorAP;
NParticlesNB = NBData.NParticlesAP;
SENBFluo = SDNBFluo./sqrt(NParticlesNB);
NBNC12 = NBData.nc12;
NBNC13 = NBData.nc13;
NBNC14 = NBData.nc14;

%% Plot the Nuclear fluorescence
AP = 12;

figure
hold on
%for AP=1:41
% nc 14
%     PlotHandle(1) = errorbar((NBnc14-31:length(NBFluo)-31)*0.66,NBFluo(NBnc14:end,AP),...
%         SENBFluo(NBnc14:end,AP))
%      PlotHandle(2) = errorbar((NBhalfnc14:length(NBhalfFluo))*0.66,NBhalfFluo(NBhalfnc14:end,AP),...
%          SENBhalfFluo(NBhalfnc14:end,AP))
%     title('Hb-NB-eGFP Nuclear fluorescence in nc 14')
%     xlabel('Time (min)') 
%     ylabel('Nuclear fluorescence (AU)')
%     leg = legend('1x dosage','0.5x dosage');
%     ylim([0 2000])
    
%    pause
%    end 
% nc 13
    PlotHandle(1) = errorbar(NBTime(NBNC13:length(NBFluo)),NBFluo(NBNC13:end,AP),...
        SENBFluo(NBNC13:end,AP))
    title('NB-NB Nuclear fluorescence in nc 13, 14')
    xlabel('Time (min)') 
    ylabel('Nuclear fluorescence (AU)')
    leg = legend('NB-NB');
    %ylim([0 2000])

%% Cytoplasmic fluo ( This should be done again with a better cytoplasmic mask)
AP = 12;
%hold on
%5plot(NBData.MeanCytoAPProfile{1,1}(:,AP))

NBCytoFluo = NBData.MeanCytoAPProfile{1,1};
NBCytoFluo = NBCytoFluo';

NBSECytoFluo = NBData.SECytoAPProfile{1,1};
NBSECytoFluo = NBSECytoFluo';
%% Cyto fluo along AP axis (I should revisit this with better nuclear mask)
hold on
for i=1:length(NBCytoFluo)
    errorbar(0:0.025:1,NBCytoFluo(i,:),NBSECytoFluo(i,:))
    pause
end
hold off 

%% Section for checking the cytoplasmic mask

%% No Nb control - Need to run this to calculate Kg 
% This analysis is from the vasa-eGFP(11);His-iRFP transheterozygote dataset
%'D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-10-11-vasa-eGFP-His-iRFP\CompiledNuclei.mat'
% This dataset has AP info
% The bleaching test has been done, I need to analyze another dataset with
% the same laser power, but different acquisition rate.
%Prefix = '2018-05-22-NB-JB3-MCP-mCherry-vasa-eGFP1'
%NoNBData = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-03-05-NoNB-MCP-mCherry-vasa-eGFP\CompiledNuclei.mat')
NoNBData = LoadMS2Sets('NoNB')
% Define AP bin to plot
APtoPlotNoNb=0.4;
[Dummy,APtoPlotNoNb]=min((NoNBData(1).APbinID-APtoPlotNoNb).^2);

% Integration area for nuclear fluorescence : I need to multiply this to
% the cytoplasmic fluo to get the ratio (Fluo_C/Fluo_N)
% This should be done better in the future, such that we fix our pipeline
% to save IntegrationArea field to CompiledNuclei.mat file.

%Prefix = '2018-03-05-NoNB-MCP-mCherry-vasa-eGFP';
Prefix = '2018-12-09-NoNB-vasa-eGFP';
load(['E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults',filesep,Prefix,filesep,[Prefix '_lin.mat']])
IntegrationArea = sum(sum(schnitzcells(1).Mask)); % This integration area is a cirle with a diameter of 2 microns. (Here, it's converted to pixels)

%This is for the specific dataset we have
 NbChannel=1;
 i=1;

%Which frames did we have good nuclear segmentation for?
%FramesToKeepNoNb=[62,64:68,70:73,75,85,86,88:100]; % For Prefix='2017-10-11-vasa-eGFP-His-iRFP'

%Cyto fluorescence
figure(1)
plot(NoNBData.ElapsedTime,...
    NoNBData.MeanCytoAPProfile{NbChannel}(APtoPlotNoNb,:)*IntegrationArea,'.-k')
title('Cytoplasmic fluorescence over time')
xlabel('Time (min)')
ylabel('Cytoplasmic Fluorescence (AU)')

%Nuclear fluorescence
figure(2)
plot(NoNBData.ElapsedTime,...
    NoNBData.MeanVectorAP(:,APtoPlotNoNb),'.-k')
title('Nuclear fluorescence over time')
xlabel('Time (min)')
ylabel('Nuclear Fluorescence (AU)')

%Get their ratio. I need to be careful with the area here
figure(3)
yRange=linspace(0,2);
PlotHandle=plot(NoNBData.ElapsedTime,...
    (NoNBData.MeanCytoAPProfile{NbChannel}(APtoPlotNoNb,:)*IntegrationArea)./...
    NoNBData.MeanVectorAP(:,APtoPlotNoNb)','.-k');
hold on
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc12)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc13)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc14)*...
    ones(size(yRange)),yRange,'--k');
hold off
ylim([0.5,2])
title('Ratio of Cytoplasm to Nucleus fluorescence')
xlabel('time (min)')
ylabel('Fluo_C/Fluo_N')




figure(4) % Fluo_C/Fluo_N for all AP bins (used nanmean)
yRange=linspace(0,2);
PlotHandle=plot(NoNBData.ElapsedTime,...
    (nanmean(NoNBData.MeanCytoAPProfile{NbChannel})*IntegrationArea)./...
    nanmean(NoNBData.MeanVectorAP'),'.-k');
hold on
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc12)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc13)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc14)*...
    ones(size(yRange)),yRange,'--k');
hold off
ylim([0.5,2])
xlabel('time (min)')
ylabel('Fluo_C/Fluo_N')


%Plot a histogram of the cytoplasmic-to-nuclear ratio
Ratios=(NoNBData.MeanCytoAPProfile{NbChannel}*IntegrationArea)./...
    NoNBData.MeanVectorAP';
figure(5)
[N,Bins] = hist(Ratios(:),50);
bar(Bins,N/sum(N))
xlim([0.5,2])
ylabel('frequency')
xlabel('Fluo_C/Fluo_N')
StandardFigure([],gca)

Kg=nanmean(Ratios(:));
SDKg=nanstd(Ratios(:));

[Kg,SDKg]

%% Check the AP dependence of No NB case
Timepoint = 20;
hold on
for i=1:length(NoNBData.ElapsedTime)
    errorbar(0:0.025:1,NoNBData.MeanVectorAP(i,:),...
        NoNBData.SDVectorAP(i,:)./sqrt(NoNBData.NParticlesAP(i,:)))
    pause
end

%% Subtract the nuclear fluo (NoNB) from nuclear fluo (w/ NB)

% First, check the NoNB nuclear fluo (over time)
NoNBNC12 = NoNBData.nc12;
NoNBNC13 = NoNBData.nc13;
NoNBNucFluo = NoNBData.MeanVectorAll;
NoNBSDFluo = NoNBData.SDVectorAll;
NoNBNParticles = NoNBData.NParticlesAll;
NoNBSEFluo = NoNBSDFluo./sqrt(NoNBNParticles);

errorbar(NoNBData.ElapsedTime,NoNBNucFluo,NoNBSEFluo)
title('Nuclear fluorescence (No NB) over time')
xlabel('Time (min)')
ylabel('Nuclear fluorescence (AU)')
legend('No NB')

%% Second, check the NB-eGFP nuclear fluo over time (over AP)
figure(1)
hold on
for i=1:41
    errorbar(NBData.ElapsedTime(NBNC13:end),NBFluo(NBNC13:end,i),SENBFluo(NBNC13:end,i))
    pause
end
title('Nuclear fluorescence over time')
xlabel('Time (min)')
ylabel('Nuclear fluorescence (AU)')


figure(2)
hold on
for i=NBNC13:length(NBFluo)
    errorbar(0:0.025:1,NBFluo(i,:),SENBFluo(i,:))
    pause
end

%% Third, subtracting the nuclear fluo (NB - noNB) 
% For now, I will use the interpolation since our two datasets have different
% frame rates
NewTime = NBData.ElapsedTime(NBNC13:end) - NBData.ElapsedTime(NBNC13);
NewNBFluo = NBFluo(NBNC13:end,:);
NewSENBFluo = SENBFluo(NBNC13:end,:);

% Use two different interpolation methods, pchip and interp1
PchipNoNBFluo = pchip((NoNBData.ElapsedTime(NoNBNC13:end)-NoNBData.ElapsedTime(NoNBNC13)),...
                NoNBNucFluo(NoNBNC13:end),NewTime);
InterpNoNBFluo = interp1((NoNBData.ElapsedTime(NoNBNC13:end)-NoNBData.ElapsedTime(NoNBNC13)),...
    NoNBNucFluo(NoNBNC13:end),NewTime);
            
hold on
plot((NoNBData.ElapsedTime(NoNBNC13:end)-NoNBData.ElapsedTime(NoNBNC13)),NoNBNucFluo(NoNBNC13:end),'-o')
plot(NewTime,PchipNoNBFluo,'-o')
plot(NewTime,InterpNoNBFluo,'-o')

% I think interp1 is better in case I have shorter NoNB traces

%%  Subtract the background now
% First, for all AP bins, let's subtract the InterpNoNBFluo in case it's
% not all NaNs.
for i=1:length(NBData.APbinID)
    % In case all the values in that AP bin are NaNs, I will leave it
    % without subtracting the background.
    if sum(isnan(NBFluo(NBNC13:end,i))) == length(NBFluo(NBNC13:end,i))
        NBFluo_BGsubtracted(:,i) = NBFluo(NBNC13:end,i);
    else
        NBFluo_BGsubtracted(:,i) = NBFluo(NBNC13:end,i) - InterpNoNBFluo' ;
    end
end
        
%% Cytoplasmic fluorescence revisit
% Issues :
% 1) Since the nuclear mask can be smaller than the acutal nuclei, the
% cytoplasmic fluorescence can be contaminated by the nuclear fluorescence,
% 2) Also, if there's reflection, this can also incrase the intensity

% Thus, I ran Simon's CytoFluo.m code which increase the mask cricle's
% radius by 45%, and get the inversed image.
% The result is saved in MeanCytoFluo.mat file, as CytoFluoDensity (Time x Z)
% This CytoFluoDensity should be multiplied by the IntegrationArea to
% get the ratio of Nuclear to Cytoplasmic fluorescence.

load('E:\Paul-J\LivemRNA\Data\DynamicsResults\2018-05-22-NB-JB3-MCP-mCherry-vasa-eGFP1\MeanCytoFluo.mat')
% Here, I need to be careful about which z-slices to include, since there
% could be reflection. (Later, I might need to choose the good time
% points, as well)
zslice1 = 2;
zslice2 = 5;
CytoFluodensity = nanmean(CytoFluoDensity(:,zslice1:zslice2),2); % Get the mean over z-slices

% YJK : I can still see the fluctuation, is this from the measurement noise?

CytoFluo = CytoFluodensity * IntegrationArea;

%Get Cytoplasmic to Nucleus fluo ratio (Fluo_C/ Fluo_N) using newly
%calculated CytoFluo with bigger mask.
% Note. I need to be careful with the IntegrationArea here
 
figure(4) % Fluo_C/Fluo_N for all AP bins (used nanmean)
yRange=linspace(0,2);
PlotHandle=plot(NoNBData.ElapsedTime,...
    CytoFluo./...
    nanmean(NoNBData.MeanVectorAP,2),'.-k');
hold on
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc12)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc13)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(NoNBData.ElapsedTime(NoNBData.nc14)*...
    ones(size(yRange)),yRange,'--k');
hold off
ylim([0.5,2])
xlabel('time (min)')
ylabel('Fluo_C/Fluo_N')

%Plot a histogram of the cytoplasmic-to-nuclear ratio
Ratios=CytoFluo./...
    nanmean(NoNBData.MeanVectorAP,2);
figure(5)
[N,Bins] = hist(Ratios(:),50);
bar(Bins,N/sum(N))
xlim([0.5,1])
ylabel('frequency')
xlabel('Fluo_C/Fluo_N')
StandardFigure([],gca)

Kg=nanmean(Ratios(:));
SDKg=nanstd(Ratios(:));

[Kg,SDKg]

%% Cytoplasmic fluorescence with NB
clear all
% Load the dataset
% First, look at the CompiledNuclei
%Prefix = '2017-12-31-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP'
Prefix = '2017-12-21-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP'
NBData = load(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix,'\CompiledNuclei.mat'])
NBCyto = load(['D:\Data\YangJoon\LivemRNA\Data\ProcessedData\',Prefix,'_\CytoImages.mat'])
load(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults',filesep,Prefix,filesep,[Prefix '_lin.mat']])

Prefix2 = '2017-12-29-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP-halfdosage'
NBHalfData = load(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix2,'\CompiledNuclei.mat'])
%NBHalfCyto = load(['D:\Data\YangJoon\LivemRNA\Data\ProcessedData\',Prefix2,'_\CytoImages.mat'])

OneX_nc13 = NBData.nc13;
OneX_nc14 = NBData.nc14;

HalfX_nc13 = NBHalfData.nc13;
HalfX_nc14 = NBHalfData.nc14;
%%
% Integration Area for nuclear fluorescence
% This integration area is a cirle with a diameter of 2 microns. (Here, it's converted to pixels)
IntegrationArea = 277;%sum(sum(schnitzcells(2).Mask)); 

% Define the Nanobody channel
NbChannel = 1;
% Define AP bin to plot
APtoPlotNb=0.3;
APbin = APtoPlotNb; % This is for the figure title (percentage)
[Dummy,APtoPlotNb]=min((NBData(1).APbinID-APtoPlotNb).^2);

%Calculate the (TF-GFP)_N using Hernan's calculation. I need to be careful with the Integration area here
% The equation is (TF-GFP)_N = Fluo_N - Fluo_C/Kg
% Here, I assumed that TF exists only in nuclei, but not in the cytoplasm,
% and Kg=0.7. I also ignored the autofluorescence for now.
Kg=0.7;

TFGFP_Nuclei(:,APtoPlotNb) = NBData.MeanVectorAP(:,APtoPlotNb)-...
    (NBData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)'/Kg;

TFGFPError(:,APtoPlotNb) = sqrt(NBData.SDVectorAP(:,APtoPlotNb).^2+...
    (NBData.SDCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)'.^2);

TFGFP_Half_Nuclei(:,APtoPlotNb) = NBHalfData.MeanVectorAP(:,APtoPlotNb)-...
    (NBHalfData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)'/Kg;

TFGFPHalfError(:,APtoPlotNb) = sqrt(NBHalfData.SDVectorAP(:,APtoPlotNb).^2+...
    (NBHalfData.SDCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)'.^2);


%%
%Cyto fluorescence
figure(1)
hold on
errorbar(NBData.ElapsedTime(OneX_nc14:end)-NBData.ElapsedTime(OneX_nc14),...
    NBData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,OneX_nc14:end)*IntegrationArea,...
    NBData.SDCytoAPProfile{NbChannel}(APtoPlotNb,OneX_nc14:end)*IntegrationArea,'.-k')

errorbar(NBHalfData.ElapsedTime(HalfX_nc14:end)-NBHalfData.ElapsedTime(HalfX_nc14),...
    NBHalfData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,HalfX_nc14:end)*IntegrationArea,...
    NBHalfData.SDCytoAPProfile{NbChannel}(APtoPlotNb,HalfX_nc14:end)*IntegrationArea,'.-r')

title(['Cytoplasmic fluorescence','AP = ',num2str(APbin*100),'%'])
xlabel('Time (min)')
ylabel('Cytoplasmic fluorescence (AU)')
legend('1x','0.5x')
hold off

%Nuclear fluorescence
figure(2)
hold on
errorbar(NBData.ElapsedTime(OneX_nc14:end)-NBData.ElapsedTime(OneX_nc14),...
    NBData.MeanVectorAP(OneX_nc14:end,APtoPlotNb),...
    NBData.SDVectorAP(OneX_nc14:end,APtoPlotNb),'.-k')

errorbar(NBHalfData.ElapsedTime(HalfX_nc14:end)-NBHalfData.ElapsedTime(HalfX_nc14),...
    NBHalfData.MeanVectorAP(HalfX_nc14:end,APtoPlotNb),...
    NBHalfData.SDVectorAP(HalfX_nc14:end,APtoPlotNb),'.-r')

title(['Nuclear fluorescence','AP = ',num2str(APbin*100),'%'])
xlabel('Time (min)')
ylabel('Cytoplasmic fluorescence (AU)')
legend('1x','0.5x')

%Calculate the (TF-GFP)_N using Hernan's calculation. I need to be careful with the area here
% The equation is (TF-GFP)_N = Fluo_N - Fluo_C/Kg
% Here, I assumed that TF exists only in nuclei, but not in the cytoplasm,
% and Kg=0.7. I also ignored the autofluorescence for now.
Kg=0.7;

figure(3)
%yRange=linspace(0,2);
hold on
errorbar(NBData.ElapsedTime(OneX_nc14:end)-NBData.ElapsedTime(OneX_nc14),...
    TFGFP_Nuclei(OneX_nc14:end,APtoPlotNb),...
    TFGFPError(OneX_nc14:end,APtoPlotNb),'.-k')

errorbar(NBHalfData.ElapsedTime(HalfX_nc14:end)-NBHalfData.ElapsedTime(HalfX_nc14),...
    TFGFP_Half_Nuclei(HalfX_nc14:end,APtoPlotNb),...
    TFGFPHalfError(HalfX_nc14:end,APtoPlotNb),'.-r')

title(['(TF-GFP)_{Nuc}','AP = ',num2str(APbin*100),'%'])
xlabel('time (min)')
ylabel('(TF-GFP)_N')
legend('1x','0.5x')
% PlotHandle=plot(NBData.ElapsedTime,...
%     (NBData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)./...
%     NBData.MeanVectorAP(:,APtoPlotNb)','.-k');
% hold on
% PlotHandle(end+1)=plot(NBData.ElapsedTime(NoNBData.nc12)*...
%     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(NBData.ElapsedTime(NoNBData.nc13)*...
%     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(NBData.ElapsedTime(NBData.nc14)*...
%     ones(size(yRange)),yRange,'--k');
% hold off
% ylim([0.5,2])


%% Plot the (TF-GFP)_Nuc for all AP bins
% First, I need to calculate the (TF-GFP)_N 
for i=1:41
    APtoPlotNb = i;
    TFGFP_Nuclei(:,APtoPlotNb) = NBData.MeanVectorAP(:,APtoPlotNb)-...
        (NBData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)'/Kg;
end

hold on
for j=1:length(NBData.ElapsedTime)
    plot(0:0.025:1, TFGFP_Nuclei(j,:))
    pause
end
%% Cytoplasmic fluorescence with NB (with a better nuclear mask, generated by Simon's code - still working...)
Prefix = '2017-12-31-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP'
NewNBCyto = load(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix,'\MeanCytoFluo.mat'])

NbChannel = 1;
% Define AP bin to plot
APtoPlotNb=0.30;
[Dummy,APtoPlotNb]=min((NBData(1).APbinID-APtoPlotNb).^2);

% Average the cytoplasmic fluorescence over z-slices, with a better mask
NewCytoFluo = nanmean(NewNBCyto.CytoFluoDensity(:,3:20),2);

%Cyto fluorescence
figure(1)
plot(NBData.ElapsedTime,...
    NewCytoFluo*IntegrationArea,'.-k')

%Nuclear fluorescence
figure(2)
plot(NBData.ElapsedTime,...
    NBData.MeanVectorAP(:,APtoPlotNb),'.-k')

figure(3)
%yRange=linspace(0,2);
TFGFP_Nuclei(:,APtoPlotNb) = NBData.MeanVectorAP(:,APtoPlotNb)-...
    (NewCytoFluo(APtoPlotNb,:)*IntegrationArea)'/Kg;

plot(NBData.ElapsedTime,...
    TFGFP_Nuclei(:,APtoPlotNb),'.-k')
% PlotHandle=plot(NBData.ElapsedTime,...
%     (NBData.MeanCytoAPProfile{NbChannel}(APtoPlotNb,:)*IntegrationArea)./...
%     NBData.MeanVectorAP(:,APtoPlotNb)','.-k');
% hold on
% PlotHandle(end+1)=plot(NBData.ElapsedTime(NoNBData.nc12)*...
%     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(NBData.ElapsedTime(NoNBData.nc13)*...
%     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(NBData.ElapsedTime(NBData.nc14)*...
%     ones(size(yRange)),yRange,'--k');
% hold off
% ylim([0.5,2])
xlabel('time (min)')
ylabel('(TF-GFP)_N')

%% Use Jacques' method for calculating the TF-GFP_Nuc
% It's basically subtracting the nuclear fluorescence of the very posterior
% nuclei from all nuclear fluorescence values, assuming that there's no Hb
% in those posterior nuclei, which is a pretty fair assumption.

% The advanatage of this approach is that we can subtract the
% autofluorescence as well as free eGFP.
% This approach assumes that the nucleus to cytoplasm free eGFP is in
% equilibrium (or quasi-equilibrium).
% Here, as a beginning, let's assume that the most posterior nuclei has no
% (or negligible) Hb protein. Then, think this value as a background fluo
% (which is combination of free eGFP + autofluorescence)
% This posterior AP = 23;

% Load the datasets
NBData = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-12-21-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP\CompiledNuclei.mat')
NBDatahalf = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-12-29-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP-halfdosage\CompiledNuclei.mat')

NBFluo = NBData.MeanVectorAP;
NBhalfFluo = NBDatahalf.MeanVectorAP;

SDNBFluo = NBData.SDVectorAP;
SDNBhalfFluo = NBDatahalf.SDVectorAP;

NParticlesNB = NBData.NParticlesAP;
NParticlesNBhalf = NBDatahalf.NParticlesAP;

SENBFluo = SDNBFluo./sqrt(NParticlesNB);
SENBhalfFluo = SDNBhalfFluo./sqrt(NParticlesNBhalf);

NBNC13 = NBData.nc13;
NBNC14 = NBData.nc14;

NBhalfnc13 = NBDatahalf.nc13;
NBhalfnc14 = NBDatahalf.nc14;

% First, let's subtract the background from the 1x dosage (eGFP)
NBFluo(isnan(NBFluo))=0; % Make Nans to 0s.
for i=1:length(NBFluo) % time points
    for j=1:size(NBFluo,2) % AP bins
        NBFluo_sub(i,j) = NBFluo(i,j) - NBFluo(i,22);
        SDNBFluo_sub(i,j) = sqrt(SDNBFluo(i,j).^2 + SDNBFluo(i,22).^2);
    end
end

% Second, 0.5x dosage(eGFP)
NBhalfFluo(isnan(NBhalfFluo))=0; % Make Nans to 0s.
for i=1:length(NBhalfFluo) % time points
    for j=1:size(NBhalfFluo,2) % AP bins
        NBhalfFluo_sub(i,j) = NBhalfFluo(i,j) - NBhalfFluo(i,22);
        SDNBhalfFluo_sub(i,j) = sqrt(SDNBhalfFluo(i,j).^2 + SDNBhalfFluo(i,22).^2);
    end
end

%% Plot to check the background subtracted fluo (which should be
% proportional to the true Hb-NB- eGFP complex concentration)
AP = 10;
hold on
    PlotHandle(1) = errorbar((NBNC14:length(NBFluo))*0.66,NBFluo_sub(NBNC14:end,AP),...
        SDNBFluo_sub(NBNC14:end,AP))
    PlotHandle(2) = errorbar((NBhalfnc14-34:length(NBhalfFluo)-34)*0.66,NBhalfFluo_sub(NBhalfnc14:end,AP),...
         SDNBhalfFluo_sub(NBhalfnc14:end,AP))
    title(['Hb-NB-eGFP Nuclear fluorescence in nc 14','@ AP = ',num2str((AP-1)*2.5),'%'])
    xlabel('Time (min)') 
    ylabel('Nuclear fluorescence (AU)')
    leg = legend('1x dosage','0.5x dosage');
    ylim([0 2000])

end