function BcdNanobody


%% Load the data sets

DataGreen=LoadMS2Sets('Bcd-Nb+GFP');
DataRed=LoadMS2Sets('Bcd-Nb+mCherry');
DataGFP=LoadMS2Sets('Bcd-GFP');
DataTitration=LoadMS2Sets('Bcd-Nb+GFP-TitrationControl');
%Bcd-GFP and Bcd-Nb+mCherry
DataGFPmCherry=LoadMS2Sets('Bcd-GFP+Bcd-Nb+mCherry');
%Autofluorescence data
DataAutoFluo=LoadMS2Sets('GFPAutoFluo');
%No nanobody
DataNoNb=LoadMS2Sets('GFP-NoNbControl');

%Bcd=GFP data sets taken with 63x
i63x=[1,2,3,7];

%Get the nuclear mask used for the protein measurements
Prefix=DataGFPmCherry(1).SetName(10:end-1);
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
IntegrationArea=sum(sum(schnitzcells(1).Mask));
clear schnitzcells



%What AP position do we want to plot?
APtoPlot=0.3;
[Dummy,APBinToPlot]=min((DataGreen(1).APbinID-APtoPlot).^2);

%Nuclear length vs. time
BrandtData=xlsread('Brandt2006-Fig1G.xlsx');


%Simulation of delay due to fluorescent protein maturation
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig1\SupportFiles\DataForDelayPlots.mat')

%Hunchback immunostaining comparison
load('ImmunoDatForHernan.mat')
ImmunoDataForHernan=ImmunoDatForHernan;

%Ftz protein traces
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig3\SupportFiles\Fig2BFtzProteinNuc161.mat')
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig3\SupportFiles\Fig2CFtzProteinNuc138Fit.mat')
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig3\SupportFiles\Fig2CFtzProteinNuc138.mat')

%Ftz MS2+protein traces
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig4\SupportFiles\Fig3CFtzProteinAndmRNANuc147.mat')
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig4\SupportFiles\Fig3DandEFtzProteinAndmRNANuc52.mat')
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig4\SupportFiles\Fig3DandEFtzProteinAndmRNANuc52Err.mat')
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig4\SupportFiles\Fig3CFtzProteinAndmRNANuc147Error.mat')


%Kruppel + eve traces
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig5\SupportFiles\Fig4DataForEve2Kr.mat')
load('C:\Users\herna\Documents\Dropbox\My Papers\NanobodiesMethodsPaper\Figures\Figures\Fig5\SupportFiles\Fig4DataForEve2KrErr.mat')



%% No Nb control - Need to run this to calculate Kg

%I'm going to have to ask for a different AP bin here, because this one
%data set was not taken at the anterior end as all other data sets! This
%shouldn't make a huge difference.
APtoPlotNoNb=0.4;
[Dummy,APtoPlotNoNb]=min((DataGreen(1).APbinID-APtoPlotNoNb).^2);

%This is for the specific data set we have
NbChannel=1;
i=1;

%Which frames did we have good nuclear segmentation for?
FramesToKeepNoNb=[11,21:28,37:length(DataNoNb(i).ElapsedTime)];


%Cyto fluorescence
figure(1)
plot(DataNoNb(i).ElapsedTime,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(APtoPlotNoNb,:),'.-k')

%Nuclear fluorescence
figure(2)
plot(DataNoNb(i).ElapsedTime,...
    DataNoNb(i).MeanVectorAP(:,APtoPlotNoNb),'.-k')


%Get their ratio. I need to be careful with the area here
figure(3)
yRange=linspace(0,2);
PlotHandle=plot(DataNoNb(i).ElapsedTime,...
    (DataNoNb(i).MeanCytoAPProfile{NbChannel}(APtoPlotNoNb,:)*sum(sum(DataNoNb(i).IntegrationArea)))./...
    DataNoNb(i).MeanVectorAP(:,APtoPlotNoNb)','.-k');
hold on
PlotHandle(end+1)=plot(DataNoNb(i).ElapsedTime(DataNoNb(i).nc12)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(DataNoNb(i).ElapsedTime(DataNoNb(i).nc13)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)*...
    ones(size(yRange)),yRange,'--k');
hold off
ylim([0.7,1.1])
xlabel('time (min)')
ylabel('Fluo_C/Fluo_N')




figure(4)
yRange=linspace(0,2);
PlotHandle=plot(DataNoNb(i).ElapsedTime,...
    (nanmean(DataNoNb(i).MeanCytoAPProfile{NbChannel})*sum(sum(DataNoNb(i).IntegrationArea)))./...
    nanmean(DataNoNb(i).MeanVectorAP'),'.-k');
hold on
PlotHandle(end+1)=plot(DataNoNb(i).ElapsedTime(DataNoNb(i).nc12)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(DataNoNb(i).ElapsedTime(DataNoNb(i).nc13)*...
    ones(size(yRange)),yRange,'--k');
PlotHandle(end+1)=plot(DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)*...
    ones(size(yRange)),yRange,'--k');
hold off
ylim([0.6,0.85])
xlabel('time (min)')
ylabel('Fluo_C/Fluo_N')


%Plot a histogram of the cytoplasmic-to-nuclear ratio
Ratios=(DataNoNb(i).MeanCytoAPProfile{NbChannel}*sum(sum(DataNoNb(i).IntegrationArea)))./...
    DataNoNb(i).MeanVectorAP';
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



%Plot the nuclear fluorescence as a function of AP for different time
%points
figure(6)
clf
hold all
ylim([0,400])
for j=1:length(DataNoNb(i).ElapsedTime)
    plot(DataNoNb(i).APbinID,...
        DataNoNb(i).MeanVectorAP(j,:),'-')
    j
    %pause
end
hold off




i=1;
NbChannel=1;
MCPChannel=2;


figure(8)
% %10 min into nc13
% FrameToPlot=27;
% PlotHandle=errorbar(DataNoNb(i).APbinID*100,...
%     DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
%     DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-g');
% DataNoNb(i).ElapsedTime(FrameToPlot)-...
%     DataNoNb(i).ElapsedTime(DataNoNb(i).nc13)
% PlotHandle(end).CapSize = 0;

%10 min into nc14
FrameToPlot=44;
PlotHandle=errorbar(DataNoNb(i).APbinID*100,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
    DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-k');
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
PlotHandle(end).CapSize = 0;
hold on
%20 min into nc14
FrameToPlot=56;
PlotHandle(end+1)=errorbar(DataNoNb(i).APbinID*100,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
    DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-g');
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
PlotHandle(end).CapSize = 0;
%40 min into nc14
FrameToPlot=75;
PlotHandle(end+1)=errorbar(DataNoNb(i).APbinID*100,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
    DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-r');
PlotHandle(end).CapSize = 0;
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
hold off
ylim([0,3.2])
xlim([0.33,0.57]*100)
set(gca,'YTick',[0:1:3])
xlabel('position along the anterior-posterior axis (%)')
ylabel('cytoplasmic fluorescence (au)')
legend('10','20','40','Location','SouthEast')
StandardFigure(PlotHandle,gca)



FrameToPlot=43;
figure(7)
PlotHandle=errorbar(DataNoNb(i).APbinID,...
    DataNoNb(i).MeanVectorAP(FrameToPlot,:),...
    DataNoNb(i).SDVectorAP(FrameToPlot,:),'.-k');
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
hold on
FrameToPlot=53;
PlotHandle(end+1)=errorbar(DataNoNb(i).APbinID,...
    DataNoNb(i).MeanVectorAP(FrameToPlot,:),...
    DataNoNb(i).SDVectorAP(FrameToPlot,:),'.-r');
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
hold off
%ylim([0,400])
xlim([0.1,0.6])
xlabel('position along the anterior-posterior axis (%)')
ylabel('nuclear fluorescence')


%% Bcd-GFP


%Plot all of them starting at nc13
figure(1)
clf
hold all
LegendString={};
for i=1:length(DataGFP)
    errorbar(DataGFP(i).ElapsedTime-DataGFP(i).ElapsedTime(DataGFP(i).nc13),...
        DataGFP(i).MeanVectorAP(:,APBinToPlot),...
        DataGFP(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(:,APBinToPlot)),'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')


%Normalized levels in nc13
figure(2)
clf
hold all
LegendString={};
for i=1:length(i63x)
    FrameRange=DataGFP(i).nc13:DataGFP(i).nc14;
    errorbar(DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc13),...
        DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)),...
        DataGFP(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(FrameRange,APBinToPlot))/...
        max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)),'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')



%Normalized levels in nc14
figure(3)
clf
hold all
LegendString={};
for i=1:length(DataGFP)
    FrameRange=DataGFP(i).nc14:length(DataGFP(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    errorbar(DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc14),...
        DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,...
        DataGFP(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(FrameRange,APBinToPlot))/...
        MaxValue,'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence/fluorescence_{max}')
ylim([0,1.1])



%Check out the Brandt data and see whether it makes any difference for our
%measurements of concentration

%Find the point where nc starts
[MinValue,BrandtIndex]=min(BrandtData(:,1).^2);
BrandtDataNC14=BrandtData(BrandtIndex:end,:);

%Normalized levels in nc14 divided by the nuclear length as reported by
%Brandt.
figure(4)
clf
hold all
LegendString={};
for i=1:length(DataGFP)-1
    %Get the frame range and the maximum value in nc14
    FrameRange=DataGFP(i).nc14:length(DataGFP(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    %Interpolate the nuclear sizes
    BrandtInterpolate=pchip(BrandtDataNC14(:,1),BrandtDataNC14(:,2),...
        DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc14));
    %Turn the first data point into NaN to avoid trouble
    BrandtInterpolate(1)=nan;

    MaxValue=MaxValue./BrandtInterpolate;
    
    
    errorbar(DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc14),...
        DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)./...
        MaxValue',...
        DataGFP(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(FrameRange,APBinToPlot))./...
        MaxValue','.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('total fluorescence')




%Normalized levels in nc14 divided multiplied the nuclear length as reported by
%Brandt.
figure(5)
clf
hold all
LegendString={};
for i=i63x
    %Get the frame range and the maximum value in nc14
    FrameRange=DataGFP(i).nc14:length(DataGFP(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot));
        
    errorbar(DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc14),...
        DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)./...
        MaxValue',...
        DataGFP(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(FrameRange,APBinToPlot))./...
        MaxValue','.-')    
    LegendString=[LegendString,num2str(i)];
end

for i=[5,6]

    %Get the frame range and the maximum value in nc14
    FrameRange=DataGFP(i).nc14:length(DataGFP(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot));

    %Interpolate the nuclear sizes
    BrandtInterpolate=pchip(BrandtDataNC14(:,1),BrandtDataNC14(:,2),...
        DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc14));
    %Turn the first data point into NaN to avoid trouble
    BrandtInterpolate(1)=nan;

    MaxValue=MaxValue.*BrandtInterpolate/6;


    errorbar(DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(DataGFP(i).nc14),...
        DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)./...
        MaxValue',...
        DataGFP(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(FrameRange,APBinToPlot))./...
        MaxValue','.-')    
    LegendString=[LegendString,num2str(i)];
end
    
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence concentration')



%Calculate an average Bcd trace in each nc
[TimeNC13GFP,TimeNC14GFP,...
    MeanVectorAPNC13InterpGFP,MeanVectorAPNC14InterpGFP,...
    MeanNC13GFP,SDNC13GFP,SENC13GFP,...
    MeanNC14GFP,SDNC14GFP,SENC14GFP]=...
    AverageInputProteinNoHistone(DataGFP(i63x),[1,1,0,0],[0,1,1,1]);


%This plot allows us to check whether the alignment was done correctly
%in the interpolation
figure(4)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpGFP)
    plot(TimeNC14GFP,MeanVectorAPNC14InterpGFP{i}(:,APBinToPlot),'.-')
end
hold off


figure(5)
clf
hold all
for i=1:length(MeanVectorAPNC13InterpGFP)
    plot(TimeNC13GFP,MeanVectorAPNC13InterpGFP{i}(:,APBinToPlot),'.-')
end
hold off


figure(6)
PlotHandle=errorbar(TimeNC14GFP,MeanNC14GFP(:,APBinToPlot),SENC14GFP(:,APBinToPlot),'.-r')
xlabel('time into nc14 (min)')
ylabel('normalized mCherry fluorescence')
xlim([0,40])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)


figure(7)
PlotHandle=errorbar(TimeNC13GFP,MeanNC13GFP(:,APBinToPlot),SENC13GFP(:,APBinToPlot),'.-r')
xlabel('time into nc13 (min)')
ylabel('normalized mCherry fluorescence')
xlim([0,15])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)

%Check the larger pinhole data
i=4;
figure(8)
clf
hold all
for j=1:length(DataGFP(i).APbinID)
    plot(DataGFP(i).ElapsedTime,DataGFP(i).MeanVectorAP(:,j),'-')
end
hold off



%% Background in Bcd-Nb-GFP

%I'm going to try to be careful about the frames where we got a nice
%nuclear segmentation.
clear FramesToKeep
FramesToKeep{1}=[1,2,8:21,29:length(DataGreen(1).ElapsedTime)];
FramesToKeep{2}=[1:9,15:length(DataGreen(2).ElapsedTime)];
FramesToKeep{3}=[1:3,8:21,28:length(DataGreen(3).ElapsedTime)];
FramesToKeep{4}=[14:19,24:37,43:length(DataGreen(4).ElapsedTime)];


%Compare the nuclear and cytoplasmic signals
figure(1)
for i=1:length(DataGreen)
    clf
    errorbar(DataGreen(i).ElapsedTime,...
        DataGreen(i).MeanVectorAP(:,APBinToPlot),...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot)),'.-g')    
    hold on
    plot(DataGreen(i).ElapsedTime,...
        DataGreen(i).MeanCytoAPProfile{1}(APBinToPlot,:)*...
        sum(sum(DataGreen(i).IntegrationArea)),'.-k')    
    plot(DataGreen(i).ElapsedTime(FramesToKeep{i}),...
        DataGreen(i).MeanCytoAPProfile{1}(APBinToPlot,FramesToKeep{i})*...
        sum(sum(DataGreen(i).IntegrationArea)),'ok')    
    errorbar(DataGreen(i).ElapsedTime(FramesToKeep{i}),...
        DataGreen(i).MeanVectorAP(FramesToKeep{i},APBinToPlot),...
        DataGreen(i).SDVectorAP(FramesToKeep{i},APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(FramesToKeep{i},APBinToPlot)),'og')    
    hold off
    box on
    %pause
end
xlabel('frame')
ylabel('raw fluorescence')
legend('nuclear','cytoplasmic')


%Is there an AP dependence in the cytoplasmic background? How does it
%compare to the nuclear fluorescence signal. Look at the maximum point in
%nc14.

i=1;

%Find the frame corresponding to the maximum in nc14
FrameRange=DataGreen(i).nc14:length(DataGreen(i).ElapsedTime);
[~,MaxIndex]=max(DataGreen(i).MeanVectorAP(:,APBinToPlot));
MaxFrame=FrameRange(MaxIndex);
MiddleFrame=MaxFrame+10;    %Look around 10 frames (~10 minutes later)

figure(2)
errorbar(DataGreen(i).APbinID,...
    DataGreen(i).MeanVectorAP(MaxFrame,:),...
    DataGreen(i).SDVectorAP(MaxFrame,:)./...
    sqrt(DataGreen(i).NParticlesAP(MaxFrame,:)),'.-g')    
hold on
plot(DataGreen(i).APbinID,...
        DataGreen(i).MeanCytoAPProfile{1}(:,MaxFrame)*...
        sum(sum(DataGreen(i).IntegrationArea)),'.-k')    
errorbar(DataGreen(i).APbinID,...
    DataGreen(i).MeanVectorAP(MiddleFrame,:),...
    DataGreen(i).SDVectorAP(MiddleFrame,:)./...
    sqrt(DataGreen(i).NParticlesAP(MiddleFrame,:)),'o-g')   
plot(DataGreen(i).APbinID,...
        DataGreen(i).MeanCytoAPProfile{1}(:,MiddleFrame)*...
        sum(sum(DataGreen(i).IntegrationArea)),'o-k')    
hold off
xlabel('AP position')
ylabel('fluorescence')
legend('nuclear','cytoplasmic')
set(gca,'YScale','log')

%How much is the change in cytoplasmic fluorescence along the AP axis?
min(DataGreen(i).MeanCytoAPProfile{1}(:,MaxFrame))/...
    max(DataGreen(i).MeanCytoAPProfile{1}(:,MaxFrame))
min(DataGreen(i).MeanCytoAPProfile{1}(:,MiddleFrame))/...
    max(DataGreen(i).MeanCytoAPProfile{1}(:,MiddleFrame))


%How much of a difference does it make to subtract the background
%fluorescence?

%Compare with and without background subtraction normalized by the maximum
%point
figure(3)

for i=1:length(DataGreen)
    MaxValue=max(DataGreen(i).MeanVectorAP(:,APBinToPlot));
    errorbar(DataGreen(i).ElapsedTime,...
        DataGreen(i).MeanVectorAP(:,APBinToPlot)/MaxValue,...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot))/MaxValue,'.-g')    
    hold on
    MaxValue=max(DataGreen(i).MeanVectorAP(:,APBinToPlot)-...
        DataGreen(i).MeanCytoAPProfile{1}(APBinToPlot,:)'*...
        sum(sum(DataGreen(i).IntegrationArea)));
    errorbar(DataGreen(i).ElapsedTime,...
        (DataGreen(i).MeanVectorAP(:,APBinToPlot)-...
        DataGreen(i).MeanCytoAPProfile{1}(APBinToPlot,:)'*...
        sum(sum(DataGreen(i).IntegrationArea)))/MaxValue,...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot))/MaxValue,'.-k')    
    hold off 
    legend('Raw','Cyto subtracted')
    xlabel('frame')
    ylabel('normalized fluorescence')
    drawnow
    pause
end

%Test the function to remove the nanobody background
DataGreenRescaled=NanobodyFluorescence(DataGreen,Kg,SDKg);

figure(4)
for i=1:length(DataGreen)
    errorbar(DataGreen(i).ElapsedTime,...
        (DataGreen(i).MeanVectorAP(:,APBinToPlot)-...
        DataGreen(i).MeanCytoAPProfile{1}(APBinToPlot,:)'*...
        sum(sum(DataGreen(i).IntegrationArea))),...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot)),'.-k')    
    hold on
    plot(DataGreenRescaled(i).ElapsedTime,...
        DataGreenRescaled(i).MeanVectorAP(:,APBinToPlot),'.-g')   
%     errorbar(DataGreenRescaled(i).ElapsedTime,...
%         DataGreenRescaled(i).MeanVectorAP(:,APBinToPlot),...
%         DataGreenRescaled(i).SEVectorAP(:,APBinToPlot),'.-g')    
    hold off 
    legend('Manually subtracted','Function subtracted')
    xlabel('frame')
    ylabel('normalized fluorescence')
    drawnow
    %pause
end


%Looking at autofluorescence in the nuclei and cytoplasm. This was taken
%using 5uW of the 488nm laser line.
%Bcd-Nb+GFP sets taken with this power
BcdNbGFp5uW=1;

APBinToPlotNoFluo=25;

figure(5)
errorbar(DataAutoFluo.ElapsedTime-...
    DataAutoFluo.ElapsedTime(DataAutoFluo.nc14),...
    DataAutoFluo.MeanVectorAP(:,APBinToPlotNoFluo),...
    DataAutoFluo.SDVectorAP(:,APBinToPlotNoFluo)./...
    sqrt(DataAutoFluo.NParticlesAP(:,APBinToPlotNoFluo)),'.-k')
hold on
plot(DataAutoFluo.ElapsedTime-...
    DataAutoFluo.ElapsedTime(DataAutoFluo.nc14),...
    DataAutoFluo.MeanCytoAPProfile{1}(APBinToPlotNoFluo,:)*IntegrationArea,'.-r')
errorbar(DataGreen(BcdNbGFp5uW).ElapsedTime-...
    DataGreen(BcdNbGFp5uW).ElapsedTime(DataGreen(BcdNbGFp5uW).nc14),...
    DataGreen(BcdNbGFp5uW).MeanVectorAP(:,APBinToPlot),...
    DataGreen(BcdNbGFp5uW).SDVectorAP(:,APBinToPlot)./...
    sqrt(DataGreen(BcdNbGFp5uW).NParticlesAP(:,APBinToPlot)),'.-g')
plot(DataGreen(BcdNbGFp5uW).ElapsedTime-...
    DataGreen(BcdNbGFp5uW).ElapsedTime(DataGreen(BcdNbGFp5uW).nc14),...
    DataGreen(BcdNbGFp5uW).MeanCytoAPProfile{1}(APBinToPlot,:)*IntegrationArea,...
    '.-b')

hold off
xlim([-5,40])
legend('AutoNuclei','AutoCyto','Bcd-Nb+GFP nuclear','Bcd-Nb+GFP nuclear cyto',...
    'Location','SouthEast')




%% Bcd-Nb+GFP 


%Plot all of them starting at nc14
figure(1)
clf
hold all
LegendString={};
for i=1:length(DataGreen)
    errorbar(DataGreen(i).ElapsedTime-DataGreen(i).ElapsedTime(DataGreen(i).nc14),...
        DataGreen(i).MeanVectorAP(:,APBinToPlot),...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot)),'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')


%Normalized levels in nc13
figure(2)
clf
hold all
LegendString={};
for i=1:length(DataGreen)
    if DataGreen(i).nc13
        FrameRange=DataGreen(i).nc13:DataGreen(i).nc14;
        errorbar(DataGreen(i).ElapsedTime(FrameRange)-DataGreen(i).ElapsedTime(DataGreen(i).nc13),...
            DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot)/...
            max(DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot)),...
            DataGreen(i).SDVectorAP(FrameRange,APBinToPlot)./...
            sqrt(DataGreen(i).NParticlesAP(FrameRange,APBinToPlot))/...
            max(DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot)),'.-')    
        LegendString=[LegendString,num2str(i)];
    end
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')



%Normalized levels in nc14, we're going to align the embryos according to
%their maximum point
figure(3)
clf
hold all
LegendString={};
for i=1:length(DataGreen)
    FrameRange=DataGreen(i).nc14:length(DataGreen(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    errorbar(DataGreen(i).ElapsedTime(FrameRange)-DataGreen(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,...
        DataGreen(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(FrameRange,APBinToPlot))/...
        MaxValue,'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')


%Calculate an average Bcd trace in each nc
[TimeNC13Green,TimeNC14Green,...
    MeanVectorAPNC13InterpGreen,MeanVectorAPNC14InterpGreen,...
    MeanNC13Green,SDNC13Green,SENC13Green,...
    MeanNC14Green,SDNC14Green,SENC14Green]=...
    AverageInputProteinNoHistone(DataGreen,[1,0,0,0],[0,1,1,1]);


%This plot allows us to check whether the alignment was done correctly
%in the interpolation
figure(4)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpGreen)
    plot(TimeNC14Green,MeanVectorAPNC14InterpGreen{i}(:,APBinToPlot),'.-')
end
hold off


figure(5)
clf
hold all
for i=1:length(MeanVectorAPNC13InterpGreen)
    plot(TimeNC13Green,MeanVectorAPNC13InterpGreen{i}(:,APBinToPlot),'.-')
end
hold off


figure(6)
PlotHandle=errorbar(TimeNC14Green,MeanNC14Green(:,APBinToPlot),SENC14Green(:,APBinToPlot),'.-b')
xlabel('time into nc14 (min)')
ylabel('normalized mCherry fluorescence')
xlim([0,40])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)


figure(7)
PlotHandle=errorbar(TimeNC13Green,MeanNC13Green(:,APBinToPlot),SENC13Green(:,APBinToPlot),'.-b')
xlabel('time into nc13 (min)')
ylabel('normalized mCherry fluorescence')
xlim([0,40])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)








%Account for the free GFP in the nucleus
DataGreenRescaled=NanobodyFluorescence(DataGreen,Kg,SDKg,CircleArea);

%Look at the cytoplasmic fluorescence
figure(8)
clf
hold all
LegendString={};
for i=1:length(DataGreenRescaled)
    FrameRange=DataGreenRescaled(i).nc14:length(DataGreenRescaled(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGreenRescaled(i).MeanCytoAPProfile{1}(APBinToPlot,FrameRange));
    
%     errorbar(DataGreenRescaled(i).ElapsedTime(FrameRange)-DataGreenRescaled(i).ElapsedTime(FrameRange(MaxIndex)),...
%         DataGreenRescaled(i).MeanVectorAP(FrameRange,APBinToPlot)/...
%         MaxValue,...
%         DataGreenRescaled(i).SDVectorAP(FrameRange,APBinToPlot)./...
%         sqrt(DataGreenRescaled(i).NParticlesAP(FrameRange,APBinToPlot))/...
%         MaxValue,'.-') 
    
    plot(DataGreenRescaled(i).ElapsedTime(FrameRange)-DataGreenRescaled(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataGreenRescaled(i).MeanCytoAPProfile{1}(APBinToPlot,FrameRange)/...
        MaxValue,'.-')    
    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')
xlim([0,40])


figure(9)
clf
hold all
LegendString={};
for i=1:length(DataGreenRescaled)
    FrameRange=DataGreenRescaled(i).nc14:length(DataGreenRescaled(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGreenRescaled(i).MeanVectorAP(FrameRange,APBinToPlot));
    
%     errorbar(DataGreenRescaled(i).ElapsedTime(FrameRange)-DataGreenRescaled(i).ElapsedTime(FrameRange(MaxIndex)),...
%         DataGreenRescaled(i).MeanVectorAP(FrameRange,APBinToPlot)/...
%         MaxValue,...
%         DataGreenRescaled(i).SDVectorAP(FrameRange,APBinToPlot)./...
%         sqrt(DataGreenRescaled(i).NParticlesAP(FrameRange,APBinToPlot))/...
%         MaxValue,'.-') 
    
    plot(DataGreenRescaled(i).ElapsedTime(FrameRange)-DataGreenRescaled(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataGreenRescaled(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,'.-')    
    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')

%% Bcd-Nb+mCherry

%Plot all of them starting at nc13
figure(1)
clf
hold all
LegendString={};
for i=1:length(DataRed)
    errorbar(DataRed(i).ElapsedTime-DataRed(i).ElapsedTime(DataRed(i).nc13),...
        DataRed(i).MeanVectorAP(:,APBinToPlot),...
        DataRed(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataRed(i).NParticlesAP(:,APBinToPlot)),'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')


%Normalized levels in nc13
figure(2)
clf
hold all
LegendString={};
for i=1:length(DataRed)
    FrameRange=DataRed(i).nc13:DataRed(i).nc14;
    errorbar(DataRed(i).ElapsedTime(FrameRange)-DataRed(i).ElapsedTime(DataRed(i).nc13),...
        DataRed(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        max(DataRed(i).MeanVectorAP(FrameRange,APBinToPlot)),...
        DataRed(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataRed(i).NParticlesAP(FrameRange,APBinToPlot))/...
        max(DataRed(i).MeanVectorAP(FrameRange,APBinToPlot)),'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')



%Normalized levels in nc14, we're going to align the embryos according to
%their maximum point
figure(3)
clf
hold all
LegendString={};
for i=1:length(DataRed)
    FrameRange=DataRed(i).nc14:length(DataRed(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataRed(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    errorbar(DataRed(i).ElapsedTime(FrameRange)-DataRed(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataRed(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,...
        DataRed(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataRed(i).NParticlesAP(FrameRange,APBinToPlot))/...
        MaxValue,'.-')    
    LegendString=[LegendString,num2str(i)];
end
hold off
box on
legend(LegendString)
xlabel('time (min)')
ylabel('fluorescence (au)')


%Calculate an average Bcd trace in each nc
[TimeNC13Red,TimeNC14Red,...
    MeanVectorAPNC13InterpRed,MeanVectorAPNC14InterpRed,...
    MeanNC13Red,SDNC13Red,SENC13Red,...
    MeanNC14Red,SDNC14Red,SENC14Red]=...
    AverageInputProteinNoHistone(DataRed);


%This plot allows us to check whether the alignment was done correctly
%in the interpolation
figure(4)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpRed)
    plot(TimeNC14Red,MeanVectorAPNC14InterpRed{i}(:,APBinToPlot),'.-')
end
hold off


figure(5)
clf
hold all
for i=1:length(MeanVectorAPNC13InterpRed)
    plot(TimeNC13Red,MeanVectorAPNC13InterpRed{i}(:,APBinToPlot),'.-')
end
hold off


figure(6)
PlotHandle=errorbar(TimeNC14Red,MeanNC14Red(:,APBinToPlot),SENC14Red(:,APBinToPlot),'.-r')
xlabel('time into nc14 (min)')
ylabel('normalized mCherry fluorescence')
xlim([0,40])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)


figure(7)
PlotHandle=errorbar(TimeNC13Red,MeanNC13Red(:,APBinToPlot),SENC13Red(:,APBinToPlot),'.-r')
xlabel('time into nc13 (min)')
ylabel('normalized mCherry fluorescence')
xlim([0,15])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)





%% Bcd-GFP data - 2-photon

[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath]=...
    DetermineLocalFolders;

%Bcd-GFP data
[StatusNum,StatusTxt]=xlsread([DropboxFolder,filesep,'DataStatus.xlsx'],'Bcd-GFP Data');

CompileRow=find(strcmp(StatusTxt(:,1),'AnalyzeLiveData Compile Nuclei'));
CompiledSets=find(strcmp(StatusTxt(CompileRow,:),'READY'));

clear SetNames
for i=1:length(CompiledSets)
    SetName=StatusTxt{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    DataBcd(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
    load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
    DataBcd(i).APDivision=APDivision;
    SetNames{i}=SetName;
end
%Now add the SetName information
for i=1:length(DataBcd)
    DataBcd(i).SetName=SetNames{i};
end


for j=1:size(DataBcd(1).MeanVectorAP,2)
    plot(DataBcd(1).ElapsedTime,DataBcd(1).MeanVectorAP(:,j))
    drawnow
    pause
end

%% Bcd-GFP vs Bcd-Nb+mCherry

%Compare the normalized fluorescence in nc14

i=1;

FrameRange=(DataGFPmCherry(i).nc14-5):length(DataGFPmCherry(i).ElapsedTime);

MaxGFP=max(DataGFPmCherry(i).MeanVectorAP{1}((FrameRange),APBinToPlot));
MaxmCherry=max(DataGFPmCherry(i).MeanVectorAP{2}((FrameRange),APBinToPlot));

figure(1)
PlotHandle=plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    DataGFPmCherry(i).MeanVectorAP{1}((FrameRange),APBinToPlot)/MaxGFP,'.-g');
hold on
PlotHandle(end+1)=plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    DataGFPmCherry(i).MeanVectorAP{2}((FrameRange),APBinToPlot)/MaxmCherry,'.-r');
hold off
xlabel('time (min)')
ylabel('normalized fluroescence')



figure(2)
plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    DataGFPmCherry(i).MeanVectorAP{2}((FrameRange),APBinToPlot),'.-r')
hold on
plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    DataGFPmCherry(i).MeanCytoAPProfile{2}(APBinToPlot,FrameRange)*IntegrationArea,'.-k')
plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    DataGFPmCherry(i).MeanVectorAP{2}((FrameRange),APBinToPlot)-...
    DataGFPmCherry(i).MeanCytoAPProfile{2}(APBinToPlot,FrameRange)'*IntegrationArea,'.-g')
hold off
set(gca,'YScale','log')


Kg=1;

figure(3)
MaxGFP=max(DataGFPmCherry(i).MeanVectorAP{1}((FrameRange),APBinToPlot));
MaxmCherry=...
    max(DataGFPmCherry(i).MeanVectorAP{2}((FrameRange),APBinToPlot)-...
    DataGFPmCherry(i).MeanCytoAPProfile{2}(APBinToPlot,FrameRange)'*IntegrationArea/Kg);

PlotHandle=plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    DataGFPmCherry(i).MeanVectorAP{1}((FrameRange),APBinToPlot)/MaxGFP,'.-g');
hold on
plot(DataGFPmCherry(i).ElapsedTime(FrameRange),...
    (DataGFPmCherry(i).MeanVectorAP{2}((FrameRange),APBinToPlot)-...
    DataGFPmCherry(i).MeanCytoAPProfile{2}(APBinToPlot,FrameRange)'*IntegrationArea/Kg)/MaxmCherry,'.-r')
hold off


%% Bcd-GFP vs. Bcd-Nb+GFP

%%%%%%%%%%%%%%%%%%%%%WARNING%%%%%%%%%%%%%%%%%%%%%%%%
%Run "No Nb control" first in order to calculate Kg%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Compare Bcd-GFP to Bcd-Nb+GFP by taking into account:
%(1) Autofluorescence in both of them. Also in the cyto case.
%(2) Subtract them.
%(3) Doing a ch2 minimization of the multiplicative factor to compare them
%rather than normalize to the maximum
%Note that here we will assume there is no Bicoid in the cytoplasm. From
%Gregor2007a (Fig. 3D), it's not clear there is that high of a concentration during
%interphase.

%Estimate the autofluorescence. I'm using Thomas' data for this, but we're
%working on having actual measurements. First, look at the raw data from
%Bcd-GFP and BcdNb-GFP at 10uW to see how different the absolute
%intensities are.
i10uWBcdGFP=3;
i10uWBcdNb=3:4;

figure(1)
clf
hold on
for i=i10uWBcdGFP
    errorbar(DataGFP(i).ElapsedTime-DataGFP(i).ElapsedTime(DataGFP(i).nc14),...
        DataGFP(i).MeanVectorAP(:,APBinToPlot),...
        DataGFP(i).SDVectorAP(:,APBinToPlot)./sqrt(DataGFP(i).NParticlesAP(:,APBinToPlot)),...
        'o-g')    
end
for i=i10uWBcdNb
    errorbar(DataGreen(i).ElapsedTime-DataGreen(i).ElapsedTime(DataGreen(i).nc14),...
        DataGreen(i).MeanVectorAP(:,APBinToPlot),...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot)),...
        '.-k')    
end
hold off
box on
xlim([-5,50])

%OK, these data sets seem to be more or less comparable. This gives me some
%confidence that, for each set regardless of power levels, I can grab a
%fraction of the total fluorescence in order to estimate my
%autofluorescence. This could potentially change the Kg measurement as
%well. I have to be careful about this!!
AutoFluoFraction=7/25;      %This is an estimate from Gregor2007b, Fig. 2B


%Calculate the average profiles. We'll introduce the autofluorescence here
%as well.




%Calculate an average Bcd-GFP trace in each nc. Note that the shifts were
%calculated in the cell "Bcd-GFP" above.
[TimeNC13GFP,TimeNC14GFP,...
    MeanVectorAPNC13InterpGFP,MeanVectorAPNC14InterpGFP,...
    MeanNC13GFP,SDNC13GFP,SENC13GFP,...
    MeanNC14GFP,SDNC14GFP,SENC14GFP]=...
    AverageInputProteinNoHistone(DataGFP(i63x),[0,0,0,0],[0,1,1,1],AutoFluoFraction,APBinToPlot);

%Calculate an average Bcd-Nb+GFP trace in each nc. Note that the shifts were
%calculated in the cell "Bcd-Nb+GFP" above.
[TimeNC13Green,TimeNC14Green,...
    MeanVectorAPNC13InterpGreen,MeanVectorAPNC14InterpGreen,...
    MeanNC13Green,SDNC13Green,SENC13Green,...
    MeanNC14Green,SDNC14Green,SENC14Green]=...
    AverageInputProteinNoHistone(DataGreen,[1,0,0,0],[0,1,1,1],AutoFluoFraction,APBinToPlot);

%Calculate an average rescaled Bcd-Nb+GFP trace in each nc. Note that the shifts were
%calculated in the cell "Bcd-Nb+GFP" above.
DataGreenRescaled=NanobodyFluorescence(DataGreen,Kg,SDKg);
[TimeNC13GreenRescaled,TimeNC14GreenRescaled,...
    MeanVectorAPNC13InterpGreenRescaled,MeanVectorAPNC14InterpGreenRescaled,...
    MeanNC13GreenRescaled,SDNC13GreenRescaled,SENC13GreenRescaled,...
    MeanNC14GreenRescaled,SDNC14GreenRescaled,SENC14GreenRescaled]=...
    AverageInputProteinNoHistone(DataGreenRescaled,[1,0,0],[0,0,1,0],AutoFluoFraction,APBinToPlot);

%Check the alignment of the green nanobody sets. The alignment of the
%Bcd-GFP ones has been checked above already.

%This plot allows us to check whether the alignment was done correctly
%in the interpolation
figure(1)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpGFP)
    plot(TimeNC14GFP,MeanVectorAPNC14InterpGFP{i}(:,APBinToPlot),'.-')
end
hold off


figure(2)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpGreenRescaled)
    plot(TimeNC14GreenRescaled,MeanVectorAPNC14InterpGreenRescaled{i}(:,APBinToPlot),'.-')
end
hold off



%Note that I still haven't figured out how to get a nice subtraction of the
%cytoplasmic fluorescence. I should get this to work by the time we submit.
figure(3)
PlotHandle=errorbar(TimeNC14GFP,MeanNC14GFP(:,APBinToPlot),SDNC14GFP(:,APBinToPlot),'.-g')
hold on
PlotHandle(end+1)=plot(TimeNC14Green-median(diff(TimeNC14Green)),MeanNC14Green(:,APBinToPlot),'.-b')
PlotHandle(end+1)=plot(TimeNC14GreenRescaled-median(diff(TimeNC14GreenRescaled)),...
    MeanNC14GreenRescaled(:,APBinToPlot),'o-b')
hold off
xlabel('time into nc14 (min)')
ylabel('normalized GFP fluorescence')
xlim([-0.5,25])
ylim([0,1.025])
set(gca,'YTick',[0:0.2:1])
legend('Bicoid-GFP','Nb nuclear raw','Nb nuclear - cyto','Location','SouthEast')
StandardFigure(PlotHandle,gca)


%Do a chi2 to figure out the best rescaling to compare Bcd-GFP and Bcd-Nb

%What index do we want to start from to compare both sets?
GFPStart=1;
NbStart=1;
%Note that Bcd-GFP is shorter
GFPToFit=MeanNC14GFP(GFPStart:end,APBinToPlot);
NbToFit=MeanNC14GreenRescaled(NbStart:length(GFPToFit)+NbStart-1,APBinToPlot);

figure(4)
plot(GFPToFit,'.-')
hold on
plot(NbToFit,'.-')
hold off

%Fit
xFit=fminsearch(@(x)chi2Bcd(GFPToFit,NbToFit,x),[1,0])


%Rescaled and offset the Nb data
Rescale=1.0807;
OffSet=0
figure(4)
PlotHandle=errorbar(TimeNC14GFP,MeanNC14GFP(:,APBinToPlot),SDNC14GFP(:,APBinToPlot),'.-g')
PlotHandle(end).CapSize=0;
hold on
PlotHandle(end+1)=errorbar(TimeNC14GreenRescaled-median(diff(TimeNC14GreenRescaled)),...
    MeanNC14GreenRescaled(:,APBinToPlot)*Rescale+OffSet,...
    SDNC14GreenRescaled(:,APBinToPlot)*Rescale,'.-b');
PlotHandle(end).CapSize=0;
hold off
xlabel('time into nc14 (min)')
ylabel('scaled eGFP fluorescence')
xlim([-0.5,25])
ylim([0,1.15])
set(gca,'YTick',[0:0.2:1])
legend('Bicoid-GFP','Bicoid-Nb+GFP','Location','SouthEast')
StandardFigure(PlotHandle,gca)




%% Compare the three measurements

close all

%nc13

figure(1)
PlotHandle=errorbar(TimeNC13GFP,MeanNC13GFP(:,APBinToPlot),SDNC13GFP(:,APBinToPlot),'.-b');
hold on
PlotHandle(end+1)=errorbar(TimeNC13Green,MeanNC13Green(:,APBinToPlot),SDNC13Green(:,APBinToPlot),'.-g');
PlotHandle(end+1)=errorbar(TimeNC13Red(2:end),MeanNC13Red(1:end-1,APBinToPlot),SDNC13Red(1:end-1,APBinToPlot),'.-r');
hold off
xlabel('time into nc13 (min)')
ylabel('normalized fluorescence')
xlim([0,20])
ylim([0,1.05])
legend('Bicoid-GFP','Bicoid-Nanobody+GFP','Bicoid-Nanobody+mCherry','Location','SouthEast')
StandardFigure(PlotHandle,gca)


%nc14
figure(2)
PlotHandle=errorbar(TimeNC14GFP(1:end-1),MeanNC14GFP(2:end,APBinToPlot),SDNC14GFP(2:end,APBinToPlot),'.-b');
hold on
PlotHandle(end+1)=errorbar(TimeNC14Green,MeanNC14Green(:,APBinToPlot),SDNC14Green(:,APBinToPlot),'.-g');
PlotHandle(end+1)=errorbar(TimeNC14Red(),MeanNC14Red(:,APBinToPlot),SDNC14Red(:,APBinToPlot),'.-r');
hold off
xlabel('time into nc14 (min)')
ylabel('normalized fluorescence')
xlim([0,30])
ylim([0,1.05])
legend('Bicoid-GFP','Bicoid-Nanobody+GFP','Bicoid-Nanobody+mCherry','Location','SouthEast')
StandardFigure(PlotHandle,gca)


%Combine both nuclear cycles
TimeShift=25;

figure(3)
PlotHandle=errorbar(TimeNC13GFP,MeanNC13GFP(:,APBinToPlot),SDNC13GFP(:,APBinToPlot),'.-b');
hold on
PlotHandle(end+1)=errorbar(TimeNC13Green,MeanNC13Green(:,APBinToPlot),SDNC13Green(:,APBinToPlot),'.-g');
PlotHandle(end+1)=errorbar(TimeNC13Red(2:end),MeanNC13Red(1:end-1,APBinToPlot),SDNC13Red(1:end-1,APBinToPlot),'.-r');
PlotHandle(end+1)=errorbar(TimeShift+TimeNC14GFP(1:end-1),MeanNC14GFP(2:end,APBinToPlot),SDNC14GFP(2:end,APBinToPlot),'.-b');
PlotHandle(end+1)=errorbar(TimeShift+TimeNC14Green,MeanNC14Green(:,APBinToPlot),SDNC14Green(:,APBinToPlot),'.-g');
PlotHandle(end+1)=errorbar(TimeShift+TimeNC14Red,MeanNC14Red(:,APBinToPlot),SDNC14Red(:,APBinToPlot),'.-r');
hold off
xlabel('time into nc13 (min)')
ylabel('normalized fluorescence')
xlim([0,60])
ylim([0,1.05])
legend('Bicoid-GFP','Bicoid-Nanobody+GFP','Bicoid-Nanobody+mCherry','Location','SouthEast')
set(gca,'XTick',[0,10,20,25,35,45,55])
set(gca,'XTickLabel',{'0','10','20','0','10',...
    '20','30'})
StandardFigure(PlotHandle,gca)
set(gcf,'Position',[730.3333  591.6667  664.6667  420.0000])


%Jacques has an interesting point about how it might not be a good idea to
%compare protein from different embryos. For example, we might have issues
%with nc having different lengths. This would presumably stretch the data
%for some embryos. One solution is to take data where the two proteins are
%present in the same embryo. This would be ideal.

%Plot the individual traces for nc14
%Normalized levels in nc14, we're going to align the embryos according to
%their maximum point
figure(4)
clf
hold on
PlotHandle=[];
for i=1:length(DataGFP)
    FrameRange=DataGFP(i).nc14:length(DataGFP(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    PlotHandle(end+1)=errorbar(DataGFP(i).ElapsedTime(FrameRange)-DataGFP(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataGFP(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,...
        DataGFP(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGFP(i).NParticlesAP(FrameRange,APBinToPlot))/...
        MaxValue,'b.-');
end
for i=1:length(DataGreen)
    FrameRange=DataGreen(i).nc14:length(DataGreen(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    PlotHandle(end+1)=errorbar(DataGreen(i).ElapsedTime(FrameRange)-DataGreen(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataGreen(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,...
        DataGreen(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataGreen(i).NParticlesAP(FrameRange,APBinToPlot))/...
        MaxValue,'g.-');
end
for i=1:length(DataRed)
    FrameRange=DataRed(i).nc14:length(DataRed(i).ElapsedTime);
    [MaxValue,MaxIndex]=max(DataRed(i).MeanVectorAP(FrameRange,APBinToPlot));
    
    PlotHandle(end+1)=errorbar(DataRed(i).ElapsedTime(FrameRange)-DataRed(i).ElapsedTime(FrameRange(MaxIndex)),...
        DataRed(i).MeanVectorAP(FrameRange,APBinToPlot)/...
        MaxValue,...
        DataRed(i).SDVectorAP(FrameRange,APBinToPlot)./...
        sqrt(DataRed(i).NParticlesAP(FrameRange,APBinToPlot))/...
        MaxValue,'r.-');
end
hold off
box on
xlabel('time (min)')
ylabel('fluorescence (au)')
xlim([-10,30])
ylim([0,1.05])
StandardFigure(PlotHandle,gca)



%% Looking at the cytoplasmic fluorescence in Nb+GFP

%Data set to look at
i=4;

%Integration area
IntegrationRadius=4;       %This is in pixels for 256x128 pixels @ 1.7 zoom. I should read this out from the MAT files.
Circle=logical(zeros(3*IntegrationRadius,3*IntegrationRadius));
Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1);
CircleArea=sum(sum(Circle));



figure(1)
PlotHandle=errorbar(DataGreen(i).ElapsedTime,DataGreen(i).MeanCytoAPProfile(APBinToPlot,:),...
    DataGreen(i).SECytoAPProfile(APBinToPlot,:),'.-k');
xlim([0,85])
StandardFigure(PlotHandle,gca)


figure(1)
PlotHandle=errorbar(DataGreen(i).ElapsedTime,DataGreen(i).MeanCytoAPProfile(APBinToPlot,:)*CircleArea/...
    mean(DataGreen(i).MeanCytoAPProfile(APBinToPlot,:)*CircleArea),...
    DataGreen(i).SECytoAPProfile(APBinToPlot,:)*CircleArea/...
    mean(DataGreen(i).MeanCytoAPProfile(APBinToPlot,:)*CircleArea),'.-k');
xlim([0,85])
StandardFigure(PlotHandle,gca)

figure(2)
PlotHandle=errorbar(DataGreen(i).ElapsedTime,DataGreen(i).MeanCytoAPProfile(APBinToPlot,:)*CircleArea,...
    DataGreen(i).SECytoAPProfile(APBinToPlot,:)*CircleArea,'.-k');
hold on
PlotHandle(end+1)=errorbar(DataGreen(i).ElapsedTime,DataGreen(i).MeanVectorAP(:,APBinToPlot),...
    DataGreen(i).SDVectorAP(:,APBinToPlot)./sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot)),'.-g');
hold off
xlabel('time (min)')
ylabel('concentration (au)')
legend('cytoplasm','nucleus')
xlim([0,85])
StandardFigure(PlotHandle,gca)





%% Nb titration control

%Divide the samples in DataTitration. I have Egfp15, Egfp25, and EgfpA
DataTitration15=DataTitration(...
    (~cellfun(@isempty,strfind(lower({DataTitration.SetName}),'egfp15')))|...
    (~cellfun(@isempty,strfind(lower({DataTitration.SetName}),'het15'))));
DataTitration25=DataTitration(...
    (~cellfun(@isempty,strfind(lower({DataTitration.SetName}),'egfp25')))|...
    (~cellfun(@isempty,strfind(lower({DataTitration.SetName}),'het25'))));
DataTitrationA=DataTitration(...
    (~cellfun(@isempty,strfind(lower({DataTitration.SetName}),'egfpa')))|...
    (~cellfun(@isempty,strfind(lower({DataTitration.SetName}),'heta'))));




%Look at the cytoplasmic levels of the different lines at the point of
%maximum fluorescence in nc14
figure(1)
clf
hold on
PlotHandle=[];
for i=1:length(DataTitration15)
    FrameRange=DataTitration15(i).nc14:length(DataTitration15(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitration15(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitration15(i).APbinID,DataTitration15(i).MeanCytoAPProfile{1}(:,MaxFrame),...
        DataTitration15(i).SDCytoAPProfile{1}(:,MaxFrame),'k.-');
end
for i=1:length(DataTitration25)
    FrameRange=DataTitration25(i).nc14:length(DataTitration25(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitration25(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitration25(i).APbinID,DataTitration25(i).MeanCytoAPProfile{1}(:,MaxFrame),...
        DataTitration25(i).SDCytoAPProfile{1}(:,MaxFrame),'ro-')
end
for i=1:length(DataTitrationA)
    FrameRange=DataTitrationA(i).nc14:length(DataTitrationA(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitrationA(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitrationA(i).APbinID,DataTitrationA(i).MeanCytoAPProfile{1}(:,MaxFrame),...
        DataTitrationA(i).SDCytoAPProfile{1}(:,MaxFrame),'gs-')
end
hold off
box on
xlabel('AP position')
ylabel('cyto fluorescence')
legend(PlotHandle([1,length(DataTitration15)+1,length(DataTitration25)+length(DataTitration15)+1]),...
    'eGFP15','eGFP25','eGFPA')
StandardFigure(PlotHandle,gca)


%Look at raw nuclear fluorescence dynamics
figure(2)
clf
hold on
PlotHandle=[];
for i=1:length(DataTitration15)
    FrameRange=DataTitration15(i).nc14:length(DataTitration15(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitration15(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitration15(i).ElapsedTime-DataTitration15(i).ElapsedTime(MaxFrame),...
        DataTitration15(i).MeanVectorAP(:,APBinToPlot),...
        DataTitration15(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataTitration15(i).NParticlesAP(:,APBinToPlot)),'k.-');
end
for i=1:length(DataTitration25)
    FrameRange=DataTitration25(i).nc14:length(DataTitration25(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitration25(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitration25(i).ElapsedTime-DataTitration25(i).ElapsedTime(MaxFrame),...
        DataTitration25(i).MeanVectorAP(:,APBinToPlot),...
        DataTitration25(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataTitration25(i).NParticlesAP(:,APBinToPlot)),'ro-');
end
for i=1:length(DataTitrationA)
    FrameRange=DataTitrationA(i).nc14:length(DataTitrationA(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitrationA(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitrationA(i).ElapsedTime-DataTitrationA(i).ElapsedTime(MaxFrame),...
        DataTitrationA(i).MeanVectorAP(:,APBinToPlot),...
        DataTitrationA(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataTitrationA(i).NParticlesAP(:,APBinToPlot)),'go-');
end
hold off
box on
xlabel('AP position')
ylabel('raw nuclear fluorescence')
legend(PlotHandle([1,length(DataTitration15)+1,length(DataTitration25)+length(DataTitration15)+1]),...
    'eGFP15','eGFP25','eGFPA')
xlim([-10,40])
StandardFigure(PlotHandle,gca)



%Nuclear fluorescence without the free GFP component.
figure(3)
clf
hold on
PlotHandle=[];
for i=1:length(DataTitration15)
    FrameRange=DataTitration15(i).nc14:length(DataTitration15(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitration15(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitration15(i).ElapsedTime-DataTitration15(i).ElapsedTime(MaxFrame),...
        DataTitration15(i).MeanVectorAP(:,APBinToPlot)-...
        DataTitration15(i).MeanCytoAPProfile{1}(APBinToPlot,:)'*...
        sum(sum(DataTitration15(i).IntegrationArea)),...
        DataTitration15(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataTitration15(i).NParticlesAP(:,APBinToPlot)),'k.-');
end
for i=1:length(DataTitration25)
    FrameRange=DataTitration25(i).nc14:length(DataTitration25(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitration25(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitration25(i).ElapsedTime-DataTitration25(i).ElapsedTime(MaxFrame),...
        DataTitration25(i).MeanVectorAP(:,APBinToPlot)-...
        DataTitration25(i).MeanCytoAPProfile{1}(APBinToPlot,:)'*...
        sum(sum(DataTitration25(i).IntegrationArea)),...
        DataTitration25(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataTitration25(i).NParticlesAP(:,APBinToPlot)),'ro-');
end
for i=1:length(DataTitrationA)
    FrameRange=DataTitrationA(i).nc14:length(DataTitrationA(i).ElapsedTime);
    [~,MaxIndex]=max(DataTitrationA(i).MeanVectorAP(FrameRange,APBinToPlot));
    MaxFrame=FrameRange(MaxIndex);
    PlotHandle(end+1)=errorbar(DataTitrationA(i).ElapsedTime-DataTitrationA(i).ElapsedTime(MaxFrame),...
        DataTitrationA(i).MeanVectorAP(:,APBinToPlot)-...
        DataTitrationA(i).MeanCytoAPProfile{1}(APBinToPlot,:)'*...
        sum(sum(DataTitrationA(i).IntegrationArea)),...
        DataTitrationA(i).SDVectorAP(:,APBinToPlot)./...
        sqrt(DataTitrationA(i).NParticlesAP(:,APBinToPlot)),'go-');
end
hold off
box on
xlabel('AP position')
ylabel('nuclear fluorescence')
legend(PlotHandle([1,length(DataTitration15)+1,length(DataTitration25)+length(DataTitration15)+1]),...
    'eGFP15','eGFP25','eGFPA')
xlim([-10,40])
StandardFigure(PlotHandle,gca)



%% Figure XX: Bicoid-GFP vs. Bicoid-Nb+GFP

%%%%%%%%%%%%%%%%%%%%%WARNING%%%%%%%%%%%%%%%%%%%%%%%%
%Run "No Nb control" first in order to calculate Kg%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Compare Bcd-GFP to Bcd-Nb+GFP by taking into account:
%(1) Autofluorescence in both of them. Also in the cyto case.
%(2) Subtract them.
%(3) Doing a ch2 minimization of the multiplicative factor to compare them
%rather than normalize to the maximum
%Note that here we will assume there is no Bicoid in the cytoplasm. From
%Gregor2007a (Fig. 3D), it's not clear there is that high of a concentration during
%interphase.

%Estimate the autofluorescence. I'm using Thomas' data for this, but we're
%working on having actual measurements. First, look at the raw data from
%Bcd-GFP and BcdNb-GFP at 10uW to see how different the absolute
%intensities are.
i10uWBcdGFP=3;
i10uWBcdNb=3:4;

figure(1)
clf
hold on
for i=i10uWBcdGFP
    errorbar(DataGFP(i).ElapsedTime-DataGFP(i).ElapsedTime(DataGFP(i).nc14),...
        DataGFP(i).MeanVectorAP(:,APBinToPlot),...
        DataGFP(i).SDVectorAP(:,APBinToPlot)./sqrt(DataGFP(i).NParticlesAP(:,APBinToPlot)),...
        'o-g')    
end
for i=i10uWBcdNb
    errorbar(DataGreen(i).ElapsedTime-DataGreen(i).ElapsedTime(DataGreen(i).nc14),...
        DataGreen(i).MeanVectorAP(:,APBinToPlot),...
        DataGreen(i).SDVectorAP(:,APBinToPlot)./sqrt(DataGreen(i).NParticlesAP(:,APBinToPlot)),...
        '.-k')    
end
hold off
box on
xlim([-5,50])

%OK, these data sets seem to be more or less comparable. This gives me some
%confidence that, for each set regardless of power levels, I can grab a
%fraction of the total fluorescence in order to estimate my
%autofluorescence. This could potentially change the Kg measurement as
%well. I have to be careful about this!!
AutoFluoFraction=7/25;      %This is an estimate from Gregor2007b, Fig. 2B


%Calculate the average profiles. We'll introduce the autofluorescence here
%as well.

%Calculate an average Bcd-GFP trace in each nc. Note that the shifts were
%calculated in the cell "Bcd-GFP" above.
[TimeNC13GFP,TimeNC14GFP,...
    MeanVectorAPNC13InterpGFP,MeanVectorAPNC14InterpGFP,...
    MeanNC13GFP,SDNC13GFP,SENC13GFP,...
    MeanNC14GFP,SDNC14GFP,SENC14GFP]=...
    AverageInputProteinNoHistone(DataGFP(i63x),[0,0,0,0],[0,1,1,1],AutoFluoFraction,APBinToPlot);

%Calculate an average Bcd-Nb+GFP trace in each nc. Note that the shifts were
%calculated in the cell "Bcd-Nb+GFP" above.
[TimeNC13Green,TimeNC14Green,...
    MeanVectorAPNC13InterpGreen,MeanVectorAPNC14InterpGreen,...
    MeanNC13Green,SDNC13Green,SENC13Green,...
    MeanNC14Green,SDNC14Green,SENC14Green]=...
    AverageInputProteinNoHistone(DataGreen,[1,0,0,0],[0,1,1,1],AutoFluoFraction,APBinToPlot);

%Calculate an average rescaled Bcd-Nb+GFP trace in each nc. Note that the shifts were
%calculated in the cell "Bcd-Nb+GFP" above.
DataGreenRescaled=NanobodyFluorescence(DataGreen,Kg,SDKg);
[TimeNC13GreenRescaled,TimeNC14GreenRescaled,...
    MeanVectorAPNC13InterpGreenRescaled,MeanVectorAPNC14InterpGreenRescaled,...
    MeanNC13GreenRescaled,SDNC13GreenRescaled,SENC13GreenRescaled,...
    MeanNC14GreenRescaled,SDNC14GreenRescaled,SENC14GreenRescaled]=...
    AverageInputProteinNoHistone(DataGreenRescaled,[1,0,0],[0,0,1,0],AutoFluoFraction,APBinToPlot);

%Check the alignment of the green nanobody sets. The alignment of the
%Bcd-GFP ones has been checked above already.

%This plot allows us to check whether the alignment was done correctly
%in the interpolation
figure(1)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpGFP)
    plot(TimeNC14GFP,MeanVectorAPNC14InterpGFP{i}(:,APBinToPlot),'.-')
end
hold off


figure(2)
clf
hold all
for i=1:length(MeanVectorAPNC14InterpGreenRescaled)
    plot(TimeNC14GreenRescaled,MeanVectorAPNC14InterpGreenRescaled{i}(:,APBinToPlot),'.-')
end
hold off


%Note that I still haven't figured out how to get a nice subtraction of the
%cytoplasmic fluorescence. I should get this to work by the time we submit.
figure(3)
PlotHandle=errorbar(TimeNC14GFP,MeanNC14GFP(:,APBinToPlot),SDNC14GFP(:,APBinToPlot),'.-g')
hold on
PlotHandle(end+1)=plot(TimeNC14Green-median(diff(TimeNC14Green)),MeanNC14Green(:,APBinToPlot),'.-b')
PlotHandle(end+1)=plot(TimeNC14GreenRescaled-median(diff(TimeNC14GreenRescaled)),...
    MeanNC14GreenRescaled(:,APBinToPlot),'o-b')
hold off
xlabel('time into nc14 (min)')
ylabel('normalized GFP fluorescence')
xlim([-0.5,25])
ylim([0,1.025])
set(gca,'YTick',[0:0.2:1])
legend('Bicoid-GFP','Nb nuclear raw','Nb nuclear - cyto','Location','SouthEast')
StandardFigure(PlotHandle,gca)


%Do a chi2 to figure out the best rescaling to compare Bcd-GFP and Bcd-Nb

%What index do we want to start from to compare both sets?
GFPStart=1;
NbStart=1;
%Note that Bcd-GFP is shorter
GFPToFit=MeanNC14GFP(GFPStart:end,APBinToPlot);
NbToFit=MeanNC14GreenRescaled(NbStart:length(GFPToFit)+NbStart-1,APBinToPlot);

figure(4)
plot(GFPToFit,'.-')
hold on
plot(NbToFit,'.-')
hold off


%Fit
xFit=fminsearch(@(x)chi2Bcd(GFPToFit,NbToFit,x),[1])

%Rescaled and offset the Nb data
Rescale=1.05;
OffSet=0;
figure(5)
PlotHandle=errorbar(TimeNC14GFP,MeanNC14GFP(:,APBinToPlot),SDNC14GFP(:,APBinToPlot),'.-g')
PlotHandle(end).CapSize=0;
hold on
PlotHandle(end+1)=errorbar(TimeNC14GreenRescaled-median(diff(TimeNC14GreenRescaled)),...
    MeanNC14GreenRescaled(:,APBinToPlot)*Rescale+OffSet,...
    SDNC14GreenRescaled(:,APBinToPlot)*Rescale,'.-b');
PlotHandle(end).CapSize=0;
hold off
xlabel('time into nc14 (min)')
ylabel('scaled eGFP fluorescence')
xlim([-0.5,25])
ylim([0,1.15])
set(gca,'YTick',[0:0.2:1])
legend('Bicoid-eGFP','LlamaTag-Bicoid+eGFP','Location','SouthEast')
StandardFigure(PlotHandle,gca)



%Find out which data sets have similar AP position ranges. I'll use this to
%draw a snapshot of Bicoid with GFP and Nb
for i=1:length(DataGFP)
    MinAPBin=min(find(sum(DataGFP(i).APFilter)>30));
    MaxAPBin=max(find(sum(DataGFP(i).APFilter)>30));
    GFP(i,:)=[i,DataGFP(i).APbinID(MinAPBin),DataGFP(i).APbinID(MaxAPBin)]
end

for i=1:length(DataGreen)
    MinAPBin=min(find(sum(DataGreen(i).APFilter)>30));
    MaxAPBin=max(find(sum(DataGreen(i).APFilter)>30));
    Green(i,:)=[i,DataGreen(i).APbinID(MinAPBin),DataGreen(i).APbinID(MaxAPBin)]
end



%Hunchback immuno comparison

%Normalize the fluorescence values by their means
HbMean=ImmunoDataForHernan.HbImmunoIntensity(:,1)/mean(ImmunoDataForHernan.HbImmunoIntensity(:,1));
NbMean=ImmunoDataForHernan.NbIntensity(:,1)/mean(ImmunoDataForHernan.NbIntensity(:,1));
HbSD=ImmunoDataForHernan.HbImmunoIntensity(:,3)/mean(ImmunoDataForHernan.HbImmunoIntensity(:,1))./...
    sqrt(ImmunoDataForHernan.HbImmunoIntensity(:,4));
NbSD=ImmunoDataForHernan.NbIntensity(:,3)/mean(ImmunoDataForHernan.NbIntensity(:,1))./...
    sqrt(ImmunoDataForHernan.NbIntensity(:,4));

%Do a linear fit
X = [ones(size(HbMean)) HbMean];
b = regress(NbMean,X)    % Removes NaN data



figure(5)
xRange=linspace(0,2);
PlotHandle=plot(HbMean,NbMean,'.k');
hold on
PlotHandle(end+1)=plot(xRange,b(1)+b(2)*xRange,'-r');
hold off
xlim([0,1.75])
ylim([0,1.75])
xlabel('normalized anti-Hunchback fluorescence')
ylabel('normalized LlamaTag-Hunchback fluorescence')
axis square
StandardFigure(PlotHandle,gca)



%I don't think this looks as good (or informative) with the error bars
figure(6)
PlotHandle=errorbarxyHG(HbMean,NbMean,HbSD,NbSD,[],[],'ok','k','MarkerFaceColor','k');
hold on
PlotHandle(end+1)=plot(xRange,b(1)+b(2)*xRange,'-r');
hold off
xlim([0,1.75])
ylim([0,1.75])
xlabel('normalized anti-Hunchback fluorescence')
ylabel('normalized LlamaTag-Hunchback fluorescence')
axis square
set(gca,'YTick',[0:0.5:2])
set(gca,'XTick',[0:0.5:2])
StandardFigure(PlotHandle,gca)

%Raw data
figure(7)
plot(ImmunoDataForHernan.HbImmunoIntensity(:,1),...
    ImmunoDataForHernan.NbIntensity(:,1),'.k')
axis square

%% Figure 1: Delay due to FP maturation

figure(1)
PlotHandle=plot(DataForDelayPlots.time,DataForDelayPlots.TotalProtein,'-k');
hold on
PlotHandle(end+1)=plot(DataForDelayPlots.time,DataForDelayPlots.Fluorescence,'-g');
hold off
xlabel('time (min)')
ylabel('number of proteins')
ylim([0,45])
xlim([0,60])
legend('total protein','fluorescent protein')
StandardFigure(PlotHandle,gca)

rr=exp(-((DataForDelayPlots.time-5)/1.5).^2)*20;
figure(2)
PlotHandle=plot(DataForDelayPlots.time,DataForDelayPlots.TotalProtein,'-k');
hold on
PlotHandle(end+1)=plot(DataForDelayPlots.time,DataForDelayPlots.Fluorescence,'-g');
PlotHandle(end+1)=plot(DataForDelayPlots.time,rr,'-r');
hold off
xlabel('time (min)')
ylabel('number of proteins')
ylim([0,45])
xlim([0,60])
legend('total protein','fluorescent protein')
StandardFigure(PlotHandle,gca)





%% Figure 2: Ftz protein translational bursts

figure(1)
PlotHandle=errorbar(Fig2BFtzProteinNuc161.TimeintoCc14inMin,Fig2BFtzProteinNuc161.ProteinCorrected/100,...
    Fig2BFtzProteinNuc161.ProteinSEM/100,'.-g')
PlotHandle(end).CapSize=0;
xlabel('time (min)')
ylabel('protein concentration (au)')
xlim([0,60])
ylim([0,1600]/100)
StandardFigure(PlotHandle,gca)

figure(2)
PlotHandle=errorbar(Fig2CFtzProteinNuc138.TimeintoCc14inMin,Fig2CFtzProteinNuc138.ProteinCorrected/100,...
    Fig2CFtzProteinNuc138.ProteinSEM/100,'.-g')
PlotHandle(end).CapSize=0;
xlabel('time (min)')
ylabel('protein concentration (au)')
xlim([0,60])
ylim([0,1600]/100)
StandardFigure(PlotHandle,gca)

figure(3)
%I'm redoing the fit in order to extend the dynamic ranges
Fit=regress(log(Fig2CFtzProteinNuc138Fit.FitY/100),[ones(size(Fig2CFtzProteinNuc138Fit.FitX)), Fig2CFtzProteinNuc138Fit.FitX]);
xRange=linspace(0,25);
PlotHandle=errorbar(Fig2CFtzProteinNuc138Fit.Time,Fig2CFtzProteinNuc138Fit.ProtCor/100,...
    Fig2CFtzProteinNuc138Fit.ProtCorSEM/100,'.g');
PlotHandle(end).CapSize=0;
hold on
PlotHandle(end+1)=plot(xRange,exp(Fit(1)+xRange*Fit(2)),'-k');
hold off
xlabel('time (min)')
ylabel('protein concentration (au)')
set(gca,'YScale','log')
xlim([0,25])
ylim([0,1001]/100)
StandardFigure(PlotHandle,gca,'Inset')


%% Figure 3: Ftz transcription and translation

RescaleFactor=500;

figure(1)
PlotHandle=errorbar(Fig3CFtzProteinAndmRNANuc147.TimeintoCc14inMin,Fig3CFtzProteinAndmRNANuc147.MS2,...
    Fig3CFtzProteinAndmRNANuc147.MS2Err,'.-r');
PlotHandle(end).CapSize=0;
hold on
PlotHandle(end+1)=errorbar(Fig3CFtzProteinAndmRNANuc147.TimeintoCc14inMin,...
    Fig3CFtzProteinAndmRNANuc147.ProteinCorrected*RescaleFactor,...
    Fig3CFtzProteinAndmRNANuc147.ProteinSEM*RescaleFactor,'.-g');
hold off
PlotHandle(end).CapSize=0;
ylim([0,7E5])
xlim([0,50])
set(gca,'XTick',[0:10:50])
xlabel('time (min)')
ylabel('number of active PolII molecules or protein concentration (au)')
legend('MS2','LlamaTag','Location','NorthWest')
StandardFigure(PlotHandle,gca)




figure(2)
PlotHandle=plot(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    sum(Fig3DandEFtzProteinAndmRNANuc52.mRNANeighbours),'.-r');
hold on
PlotHandle(end+1)=errorbar(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    Fig3DandEFtzProteinAndmRNANuc52.ProteinCorrected*RescaleFactor,...
    Fig3DandEFtzProteinAndmRNANuc52.ProteinSEM*RescaleFactor,'.-g');
PlotHandle(end).CapSize=0;
hold off
ylim([0,1E6])
xlim([0,50])
set(gca,'XTick',[0:10:50])
xlabel('time (min)')
ylabel('number of active PolII molecules or protein concentration (au)')
legend('MS2','LlamaTag','Location','NorthEast')
StandardFigure(PlotHandle,gca)

figure(3)
PlotHandle=plot(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    sum(Fig3DandEFtzProteinAndmRNANuc52.mRNANeighbours),'.-r');
hold off
ylim([0,1E6])
xlim([0,50])
set(gca,'XTick',[0:10:50])
xlabel('time (min)')
ylabel('number of active PolII molecules (au)')
StandardFigure(PlotHandle,gca)

figure(4)
PlotHandle=plot(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    sum(Fig3DandEFtzProteinAndmRNANuc52.Nuc52mRNA),'.-r');
hold on
PlotHandle(end+1)=errorbar(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    Fig3DandEFtzProteinAndmRNANuc52.ProteinCorrected*RescaleFactor,...
    Fig3DandEFtzProteinAndmRNANuc52.ProteinSEM*RescaleFactor,'.-g');
PlotHandle(end).CapSize=0;
hold on
hold off
ylim([0,4E5])
xlim([0,50])
set(gca,'XTick',[0:10:50])
xlabel('time (min)')
ylabel('protein concentration (au)')
StandardFigure(PlotHandle,gca)




figure(5)
PlotHandle=errorbar(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    sum(Fig3DandEFtzProteinAndmRNANuc52.mRNANeighbours),...
    Fig3DandEFtzProteinAndmRNANuc52.MS2Err,'.-r');
PlotHandle(end).CapSize=0;
hold on
PlotHandle(end+1)=plot(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    Fig3DandEFtzProteinAndmRNANuc52.Nuc52mRNA,'o-r');
PlotHandle(end+1)=errorbar(Fig3DandEFtzProteinAndmRNANuc52.TimeintoCc14inMin,...
    Fig3DandEFtzProteinAndmRNANuc52.ProteinCorrected*RescaleFactor,...
    Fig3DandEFtzProteinAndmRNANuc52.ProteinSEM*RescaleFactor,'.-g');
PlotHandle(end).CapSize=0;
hold on
hold off
ylim([0,1E6])
xlim([0,50])
set(gca,'XTick',[0:10:50])
set(gca,'YTick',[0:2:10]*1E5)
xlabel('time (min)')
ylabel('protein concentration (au)')
StandardFigure(PlotHandle,gca)


%% Figure 4: Kruppel-eve input-output function

RescaleFactor=2/1000;

figure(1)
PlotHandle=errorbar(DataForEve2Kr.Time,DataForEve2Kr.ProtPlot1*RescaleFactor,...
    DataForEve2Kr.ProtPlot1SE*RescaleFactor,'.-g')
PlotHandle.CapSize=0;
hold on
PlotHandle(end+1)=errorbar(DataForEve2Kr.Time,...
    DataForEve2Kr.mRNAPlot1*RescaleFactor,DataForEve2Kr.mRNAPlot1Err*RescaleFactor,'.-r');
PlotHandle(end).CapSize=0;
hold off
xlabel('time (min)')
ylabel('fluorescence (au)')
ylim([0,10])
xlim([0,35])
set(gca,'YTick',[0:2:10])
set(gca,'XTick',[0:10:40])
StandardFigure(PlotHandle,gca)


figure(2)
PlotHandle=errorbar(DataForEve2Kr.Time,DataForEve2Kr.ProtPlot2*RescaleFactor,...
    DataForEve2Kr.ProtPlot2SE*RescaleFactor,'.-g')
PlotHandle.CapSize=0;
hold on
PlotHandle(end+1)=errorbar(DataForEve2Kr.Time,...
    DataForEve2Kr.mRNAPlot2*RescaleFactor,...
    DataForEve2Kr.mRNAPlot2Err*RescaleFactor,'.-r');
PlotHandle(end).CapSize=0;
hold off
xlabel('time (min)')
ylabel('fluorescence (au)')
ylim([0,10])
xlim([0,35])
set(gca,'YTick',[0:2:10])
set(gca,'XTick',[0:10:40])
StandardFigure(PlotHandle,gca)



%% Figure S1: cytoplasmic levels of GFP in the absence of nanobody

%Data from Jacques and Matty. There's still a little bit of concern of
%mCherry bleading into the GFP channel here. They're looking into taking
%new data using the iRFP-Histone instead of mCherry.

close all

i=1;
NbChannel=1;
MCPChannel=2;



figure(1)
%10 min into nc14
FrameToPlot=44;
PlotHandle=errorbar(DataNoNb(i).APbinID*100,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
    DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-k');
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
PlotHandle(end).CapSize = 0;
hold on
%20 min into nc14
FrameToPlot=56;
PlotHandle(end+1)=errorbar(DataNoNb(i).APbinID*100,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
    DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-g');
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
PlotHandle(end).CapSize = 0;
%40 min into nc14
FrameToPlot=75;
PlotHandle(end+1)=errorbar(DataNoNb(i).APbinID*100,...
    DataNoNb(i).MeanCytoAPProfile{NbChannel}(:,FrameToPlot),...
    DataNoNb(i).SDCytoAPProfile{NbChannel}(:,FrameToPlot),'.-r');
PlotHandle(end).CapSize = 0;
DataNoNb(i).ElapsedTime(FrameToPlot)-...
    DataNoNb(i).ElapsedTime(DataNoNb(i).nc14)
hold off
ylim([0,3.2])
xlim([0.33,0.57]*100)
set(gca,'YTick',[0:1:3])
xlabel('position along the anterior-posterior axis (%)')
ylabel('cytoplasmic fluorescence (au)')
legend('10','20','40','Location','SouthEast')
StandardFigure(PlotHandle,gca)


%Plot a histogram of the cytoplasmic-to-nuclear ratio
figure(2)
Ratios=(DataNoNb(i).MeanCytoAPProfile{NbChannel}*sum(sum(DataNoNb(i).IntegrationArea)))./...
    DataNoNb(i).MeanVectorAP';
figure(5)
[N,Bins] = hist(Ratios(:),50);
bar(Bins,N/sum(N))
xlim([0.5,1])
ylabel('frequency')
xlabel('Fluo_C/Fluo_N')
xlim([0.65,1])
StandardFigure([],gca)

Kg=nanmean(Ratios(:));
SDKg=nanstd(Ratios(:));

[Kg,SDKg]


%% Figure S2: LlamaTagRed


%Find out which data sets have similar AP position ranges. I'll use this to
%draw a snapshot of Bicoid with GFP and Nb
for i=1:length(DataGFP)
    MinAPBin=min(find(sum(DataGFP(i).APFilter)>30));
    MaxAPBin=max(find(sum(DataGFP(i).APFilter)>30));
    GFP(i,:)=[i,DataGFP(i).APbinID(MinAPBin),DataGFP(i).APbinID(MaxAPBin)];
end

for i=1:length(DataGreen)
    MinAPBin=min(find(sum(DataGreen(i).APFilter)>30));
    MaxAPBin=max(find(sum(DataGreen(i).APFilter)>30));
    Green(i,:)=[i,DataGreen(i).APbinID(MinAPBin),DataGreen(i).APbinID(MaxAPBin)];
end

for i=1:length(DataRed)
    MinAPBin=min(find(sum(DataRed(i).APFilter)>30));
    MaxAPBin=max(find(sum(DataRed(i).APFilter)>30));
    Red(i,:)=[i,DataRed(i).APbinID(MinAPBin),DataRed(i).APbinID(MaxAPBin)];
end

%Go for:
%GFP: 'Prefix ='2017-02-23-Bcd-GFP'
%Green: Prefix ='2017-02-28-BcdNbGFPA'
%Red: Prefix ='2016-12-05-BcdNbmCherryV2'

