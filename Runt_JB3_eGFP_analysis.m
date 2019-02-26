% RuntN-JB3 x vasa-eGFP1 data analysis

% Author : Yang Joon Kim
% Last Update : Dec, 2018
% Description : This code is to analyze the JB3-Runt (N-terminal fusion) x
% ubiquitous eGFP(vasa-eGFP1)


% Load the dataset
%clear all
%RuntData = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-04-29-Runt-JB3-MCP-mCherry-vasa-eGFP1-edited-1percentAPbins\CompiledNuclei.mat')
%RuntData = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-03-15-NoNB-vasa-mCherry3B-His-iRFP-transhets-nc12-14\CompiledNuclei.mat')
% Oct/2018, Now that I have homozygous embryos with Runt (or hemizygous)
RuntData = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-10-07-hbP2-r3-MS2V5-lacZ-RuntJB3-2xMCP-mCherry-vasa-eGFP1\CompiledNuclei.mat')
%% Nuclear fluorescence over AP (over time)
% Here, I am plotting the nuclear fluo over AP for all time points.
 
% Note. The background(autofluorescence) is not considered in here.
% I need a separate dataset where there was no NB, then I can get the ratio
% of free mCherry between cytoplasm and nucleus, K_mCh (equilibrium
% constant)

Time = RuntData.ElapsedTime;
MeanFluo = RuntData.MeanVectorAP;
SDFluo = RuntData.SDVectorAP;
NParticles = RuntData.NParticlesAP;
SEFluo = SDFluo./sqrt(NParticles);

%Color(gradation from blue to red) for each time point
iStart=1; %Start time
iEnd=length(Time); %End time
colormap(jet(256));
cmap=colormap;
Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
%figure(2)
hold on
%ylim([min(min(MeanFluo)) max(max(MeanFluo))])
% beginning of nc13 = 19;
% mid of nc14 = 60
if RuntData.nc12==0
    TimeVector = [RuntData.nc13,RuntData.nc14,length(RuntData.ElapsedTime)];
else
    TimeVector = [RuntData.nc12+20,RuntData.nc13+20,RuntData.nc14+20,length(RuntData.ElapsedTime)];
end
for i=TimeVector
    errorbar(0:0.025:1,MeanFluo(i,:),SEFluo(i,:),'color',Color(i-iStart+1,:))
    ylim([0 400])
    title('Nuclear Fluorescen along AP')
    xlabel('AP')
    ylabel('Nuclear Fluorescence (AU)')
    pause
end
legend('mid nc12','mid nc13','early nc14','late nc14')
standardizeFigure(gca,legend,[])

%% Nuclear Fluorescence over time (at specific AP bins)
% Here, I want to plot the Nuc. Fluo over time, inside the stripe or
% outside of the stripe region.
AP_stripe = 0.325; %0.325
AP1=AP_stripe;

[Dummy,AP_stripe]=min((RuntData.APbinID-AP_stripe).^2);

AP_nonstripe = 0.55; 
AP2=AP_nonstripe;
[Dummy,AP_nonstripe]=min((RuntData.APbinID-AP_nonstripe).^2);

figure(2)
hold on
errorbar(Time(RuntData.nc12:end)-Time(RuntData.nc12),MeanFluo(RuntData.nc12:end,AP_stripe),SEFluo(RuntData.nc12:end,AP_stripe))
errorbar(Time(RuntData.nc12:end)-Time(RuntData.nc12),MeanFluo(RuntData.nc12:end,AP_nonstripe),SEFluo(RuntData.nc12:end,AP_nonstripe))
title('Nuclear Fluorescence over Time')
xlabel('Time (min)')
ylabel('Nuclear Fluorescence (AU)')
legend(['Stripe @',num2str(AP1*100),'%'],['non-Stripe @',num2str(AP2*100),'%'])
standardizeFigure(gca,legend,[])
line([Time(RuntData.nc12)-Time(RuntData.nc12) Time(RuntData.nc12)-Time(RuntData.nc12)],[0 400],'Color','k','LineStyle','--')
line([Time(RuntData.nc13)-Time(RuntData.nc12) Time(RuntData.nc13)-Time(RuntData.nc12)],[0 400],'Color','k','LineStyle','--')
line([Time(RuntData.nc14)-Time(RuntData.nc12) Time(RuntData.nc14)-Time(RuntData.nc12)],[0 400],'Color','k','LineStyle','--')


%% Stripe region Nuclear Fluo
AP_stripe = [0.31 0.40 0.48 0.56 0.63];

hold on
for i=1:length(AP_stripe)
    AP1(i)=AP_stripe(i);
    [Dummy,AP_stripe(i)]=min((RuntData.APbinID-AP_stripe(i)).^2);
    errorbar(Time,MeanFluo(:,AP_stripe(i)),SEFluo(:,AP_stripe(i)))
    title('Nuclear Fluorescence over Time')
    xlabel('Time (min)')
    ylabel('Nuclear Fluorescence (AU)')
    ylim([1600 2400])
    pause
end
%     title('Nuclear Fluorescence over Time')
%     xlabel('Time (min)')
%     ylabel('Nuclear Fluorescence (AU)')
%     ylim([1600 2400])
    legend = [legend,['Stripe @',num2str(AP1(i)*100),'%'] ]
    
%% Interstripe region Nuclear Fluo
AP_interstripe = [0.27 0.36 0.44 0.53 0.59];

hold on
for i=1:length(AP_interstripe)
    AP2(i)=AP_interstripe(i);
    [Dummy,AP_interstripe(i)]=min((RuntData.APbinID-AP_interstripe(i)).^2);
    errorbar(Time,MeanFluo(:,AP_interstripe(i)),SEFluo(:,AP_interstripe(i)))
    pause
end
    title('Nuclear Fluorescence over Time')
    xlabel('Time (min)')
    ylabel('Nuclear Fluorescence (AU)')
    ylim([1600 2400])
    %legend = [legend,['Stripe @',num2str(AP1(i)*100),'%'] ]

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Old script for RuntN-JB10 x vasa-mCherry(3B),hets
 %% Cytoplasmic fluorescence
% 
% RuntSchnitz = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-03-15-Runt-JB10-vasa-mCh3B-His-iRFP-transhets\2018-03-15-Runt-JB10-vasa-mCh3B-His-iRFP-transhets_lin.mat')
% RuntSchnitz = RuntSchnitz.schnitzcells;
% IntegrationArea = sum(sum(RuntSchnitz(5).Mask));
% 
% CytoFluo = RuntData.MeanCytoAPProfile{1}*IntegrationArea;
% SDCytoFluo = RuntData.SDCytoAPProfile{1}*IntegrationArea;
% figure(1)
% hold on
% %ylim([min(min(MeanFluo)) max(max(MeanFluo))])
% for i=1:length(Time)
%     plot(0:0.05:1,CytoFluo(:,i),'color',Color(i-iStart+1,:))
%     pause
% end

%% Cytoplasmic fluo calculation
% Here, based on the previous plot, I will assume that the cytoplasmic
% fluorescence doesn't have AP dependence. This assumption makes sense
% since we assume the embryo as a reservoir contains many nuclei floating.

% Take the average of cytoplasmic fluo over all AP bins, then get a number
% for each time frame.
CytoFluo = nanmean(CytoFluo);
SDCytoFluo = nanmean(SDCytoFluo);
% %% No Nb control - Need to run this to calculate K_mCherry
% 
% % Load the dataset without the NB. This should change as soon as I acquire
% % a clearly noNB dataset.
% %DataNoNb = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-02-05-RuntN5-vasa-mCherry3B-His-iRFP\CompiledNuclei.mat')
% DataNoNb = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-03-15-NoNB-vasa-mCherry3B-His-iRFP-transhets-nc12-14\CompiledNuclei.mat')
% %I'm going to have to ask for a different AP bin here, because this one
% %data set was not taken at the anterior end as all other data sets! This
% %shouldn't make a huge difference.
% % Define AP bin to plot
% APtoPlotNoNb=0.4;
% [Dummy,APtoPlotNoNb]=min((DataNoNb(1).APbinID-APtoPlotNoNb).^2);
% 
% % Integration area for nuclear fluorescence : I need to multiply this to
% % the cytoplasmic fluo to get the ratio (Fluo_C/Fluo_N)
% Prefix = '2018-02-05-RuntN5-vasa-mCherry3B-His-iRFP'
% load(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults',filesep,Prefix,filesep,[Prefix '_lin.mat']])
% IntegrationArea = sum(sum(schnitzcells(5).Mask)); % This integration area is a cirle with a diameter of 2 microns. (Here, it's converted to pixels)
% 
% %This is for the specific dataset we have
%  NbChannel=1;
%  i=1;
% 
% %Which frames did we have good nuclear segmentation for?
% %FramesToKeepNoNb=[62,64:68,70:73,75,85,86,88:100]; % For Prefix='2017-10-11-vasa-eGFP-His-iRFP'
% 
% %Cyto fluorescence
% figure(1)
% plot(DataNoNb.ElapsedTime,...
%     DataNoNb.MeanCytoAPProfile{NbChannel}(APtoPlotNoNb,:),'.-k')
% title('Cytoplasmic fluorescence over time')
% xlabel('Time (min)')
% ylabel('Cytoplasmic Fluorescence (AU)')
% 
% %Nuclear fluorescence
% figure(2)
% plot(DataNoNb.ElapsedTime,...
%     DataNoNb.MeanVectorAP(:,APtoPlotNoNb),'.-k')
% title('Nuclear fluorescence over time')
% xlabel('Time (min)')
% ylabel('Nuclear Fluorescence (AU)')
% 
% %Get their ratio. I need to be careful with the area here
% figure(3)
% yRange=linspace(0,2);
% PlotHandle=plot(DataNoNb.ElapsedTime,...
%     (DataNoNb.MeanCytoAPProfile{NbChannel}(APtoPlotNoNb,:)*IntegrationArea)./...
%     DataNoNb.MeanVectorAP(:,APtoPlotNoNb)','.-k');
% hold on
% % PlotHandle(end+1)=plot(DataNoNb.ElapsedTime(DataNoNb.nc12)*...
% %     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(DataNoNb.ElapsedTime(DataNoNb.nc13)*...
%     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(DataNoNb.ElapsedTime(DataNoNb.nc14)*...
%     ones(size(yRange)),yRange,'--k');
% hold off
% ylim([0.5,2])
% title('Ratio of Cytoplasm to Nucleus fluorescence')
% xlabel('time (min)')
% ylabel('Fluo_C/Fluo_N')
% 
% 
% 
% 
% figure(4) % Fluo_C/Fluo_N for all AP bins (used nanmean)
% yRange=linspace(0,2);
% PlotHandle=plot(DataNoNb.ElapsedTime,...
%     (nanmean(DataNoNb.MeanCytoAPProfile{NbChannel})*IntegrationArea)./...
%     nanmean(DataNoNb.MeanVectorAP'),'.-k');
% hold on
% % PlotHandle(end+1)=plot(DataNoNb.ElapsedTime(DataNoNb.nc12)*...
% %     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(DataNoNb.ElapsedTime(DataNoNb.nc13)*...
%     ones(size(yRange)),yRange,'--k');
% PlotHandle(end+1)=plot(DataNoNb.ElapsedTime(DataNoNb.nc14)*...
%     ones(size(yRange)),yRange,'--k');
% hold off
% ylim([0.5,2])
% xlabel('time (min)')
% ylabel('Fluo_C/Fluo_N')
% 
% 
% %Plot a histogram of the cytoplasmic-to-nuclear ratio
% Ratios=(DataNoNb.MeanCytoAPProfile{NbChannel}*IntegrationArea)./...
%     DataNoNb.MeanVectorAP';
% figure(5)
% [N,Bins] = hist(Ratios(:),50);
% bar(Bins,N/sum(N))
% xlim([0.5,2])
% ylabel('frequency')
% xlabel('Fluo_C/Fluo_N')
% title('K_{mCh} histogram')
% StandardFigure([],gca)
% 
% Km=nanmean(Ratios(:));
% SDKm=nanstd(Ratios(:));
% 
% [Km,SDKm]
% % Note. (3/1/2018, YJK)
% % From this analysis, it seems that the K_mCherry~0.74, with SD=0.03
% % Put these plots to somewhere, like an overleaf document.
% 
% % Note that I need to subtract the autofluorescence from His-iRFP only
% % embryo, which is taken for with the red laser as well.
% 
% % For that, I need to repeat this experiment, 
% % 1) vasa-mCherry 3B hets, not fertilized with Runt-NB (should have
% % His-iRFP for nuclear marker)
% % 2) His-iRFP without the vasa-mCherry with exactly the same laser power as 1), to subtract the autofluorescence.
% % All of the detailed info for which imaging I would need is in either
% % OneNOte-Opposing gradients or my labnote.
% 
% %% Autofluorescence dataset
% % This nuclear fluo / cyto fluo should be subtracted from all of the
% % datasets
% % Here, I will start from one dataset that I took for His-iRFP only (hets)
% % Load the dataset
% AutoFluoData = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-03-19-His-iRFP-hets-autofluo\CompiledNuclei.mat')
% 
% AutoFluo = AutoFluoData.MeanVectorAP;
% SDAutoFluo = AutoFluoData.SDVectorAP;
% AutoFluo_Cyto = AutoFluoData.MeanCytoAPProfile{1};
% IntegrationArea_Auto = sum(sum(schnitzcells(5).Mask)); % This integration area is a cirle with a diameter of 2 microns. (Here, it's converted to pixels)
% AutoFluo_Cyto = AutoFluo_Cyto'.*IntegrationArea_Auto;
% 
% % Check the nuclear fluorescence
% hold on
% for j=1:size(AutoFluo,1)
%     errorbar(0:0.01:1,AutoFluo(j,:),SDAutoFluo(j,:))
%     plot(0:0.01:1,AutoFluo_Cyto(j,:))
%     pause
% end
% 
% % It seems that the autofluorescence accoutns for less than 2% of the
% % nuclear fluorescence. Thus, let's ignore this for now.
% 
% %% Subtracting the free mCherry from the Nuc Fluo (working)
% 
% %IntegrationArea = 409;%sum(sum(schnitzcells(2).Mask)); 
% NBData =  load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-03-15-Runt-JB10-vasa-mCh3B-His-iRFP-transhets\CompiledNuclei.mat')
% % Define the Nanobody channel
% NbChannel = 1;
% % Define AP bin to plot
% APtoPlotNb=0.4;
% APbin = APtoPlotNb; % This is for the figure title (percentage)
% [Dummy,APtoPlotNb]=min((NBData(1).APbinID-APtoPlotNb).^2);
% 
% % Define the Integration Area using the schnitzcells (nuclear mask)
% 
% %Calculate the (TF-GFP)_N using Hernan's calculation. I need to be careful with the Integration area here
% % The equation is (TF-GFP)_N = Fluo_N - Fluo_C/Kg
% % Here, I assumed that TF exists only in nuclei, but not in the cytoplasm,
% % and KmCh=0.74. I also ignored the autofluorescence for now.
% % From '2018-03-19-His-iRFP-hets-autofluo' dataset, nuclear autofluorescence is
% % pretty consistent over AP and time,which is around 40.
% 
% Km=0.74;
% for APtoPlotNb = 1:101
%     TF_Nuclei(:,APtoPlotNb) = NBData.MeanVectorAP(:,APtoPlotNb)-...
%         CytoFluo'/Km;
% 
%     TFError(:,APtoPlotNb) = sqrt(NBData.SDVectorAP(:,APtoPlotNb).^2+...
%         (SDCytoFluo/Km)'.^2);
%     
%     TFError(:,APtoPlotNb) = TFError(:,APtoPlotNb)./sqrt(NParticles(:,APtoPlotNb));
% end
% hold on
% for i=1:50
%     errorbar(0:0.01:1,TF_Nuclei(i,:),TFError(i,:))
%     pause
% end
% %% No NB nuclear fluo over AP
% hold on
% for i=1:length(DataNoNb.MeanVectorAP(:,1))
%     plot(0:0.01:1,DataNoNb.MeanVectorAP(i,:))
%     pause
%     
% end
%     
% %% Another approach to get TF concentration : subtracting the nuclear fluorescence from the embryo in the absence of NB
% mCh_free_Nuc = DataNoNb.MeanVectorAll(DataNoNb.nc14-3:end);
% mCh_free_Nuc_Error = DataNoNb.SDVectorAll(DataNoNb.nc14-3:end);
% 
% for i=1:length(RuntData.APbinID)
%     TF_mCh(:,i) = RuntData.MeanVectorAP(:,i) - mCh_free_Nuc';
%     TF_mCh_Error(:,i) = sqrt(RuntData.SDVectorAP(:,i).^2 + mCh_free_Nuc_Error'.^2 )./...
%                         sqrt(RuntData.NParticlesAP(:,i));
% end
%% save the analysis
% save('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\OpposingGradientAnalysis\Runt-JB10-mCherry-2018-03-15.mat',...
%     'IntegrationArea','Km','NBData','RuntSchnitz','MeanFluo','SDFluo','SEFluo','CytoFluo','SDCytoFluo','Time','TF_Nuclei','TFError')