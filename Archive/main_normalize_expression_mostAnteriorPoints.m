function main_normalize_expression_mostAnteriorPoints
%% Description 
% This script is for exploring the fold-change idea of 17.5% of APbin as a
% proxy for Runt nulls. 

%% Scheme
% 1) Take 17.5% of APbin's mean spot fluo data for each dataset.
% 2) Divide the mean spot fluo from each AP bin with 1), either for the
% whole time or for the initial slope only
% 3) First, compile all datasets with their sex
% Start with [000], [111], then move to [100]/[010]/[001]
% Caveat : Pick female embryos only at this moment.

%% (Optional) Average datasets (processing for synchronization)
% For now, I'll only focus on [000], [111], and [100]/[010]/[001] variants.

% file path to save the files
DropboxPath = 'S:/YangJoon/Dropbox';
filePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\TxnOutput_sexed';

% AverageDatasets('r1-close-female','NC',13,'savePath',filePath)
% AverageDatasets('r1-mid-female','NC',13,'savePath',filePath)
%% Load the datasets
% [000]
Data_r0 = load([filePath, filesep, 'r0-new.mat']);
% 1 Runt site
Data_r1 = load([filePath, filesep, 'r1-new-female.mat']);
Data_r1_close = load([filePath, filesep, 'r1-close-female.mat']);
Data_r1_mid = load([filePath, filesep, 'r1-mid-female.mat']);

% [111]
Data_r3 = load([filePath, filesep, 'r3-new-female.mat']);
%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
purple = [171,133,172]/255;
colorDict.purple = (4*purple - [1,1,1])/3;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
brown = [179,155,142]/255;
colorDict.brown = (2*brown - [1,1,1])/1;

colorDict.darkgreen = [126,157,144]/255;
colorDict.lightgreen = [205,214,209]/255;
thickpink = [132,27,69]/255;
colorDict.thickpink = (3*thickpink + [1,1,1]) / 4; % adding white

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.magenta; colorDict.thickpink];%;...
                %colormap(10,:); colormap(14,:); colormap(12,:); colormap(11,:)]; 
            
%% Check the MeanVectorAP for [000]
nc13_r0 = Data_r0.nc13;
nc14_r0 = Data_r0.nc14;
Time_r0 = Data_r0.ElapsedTime;
MeanFluo_r0 = Data_r0.MeanVectorAP;
SEFluo_r0 = Data_r0.SEVectorAP;

MeanFluo_individual_r0 = Data_r0.MeanVectorAP_individual;
%% Check the MeanVectorAP for [100]
nc13_r1 = Data_r1.nc13;
nc14_r1 = Data_r1.nc14;
Time_r1 = Data_r1.ElapsedTime;
MeanFluo_r1 = Data_r1.MeanVectorAP;
SEFluo_r1 = Data_r1.SEVectorAP;

MeanFluo_individual_r1 = Data_r1.MeanVectorAP_individual;

%% [010] : r1-mid
nc13_r1_mid = Data_r1_mid.nc13;
nc14_r1_mid = Data_r1_mid.nc14;
Time_r1_mid = Data_r1_mid.ElapsedTime;
MeanFluo_r1_mid = Data_r1_mid.MeanVectorAP;
SEFluo_r1_mid = Data_r1_mid.SEVectorAP;

MeanFluo_individual_r1_mid = Data_r1_mid.MeanVectorAP_individual;

%% [001] : r1-close
nc13_r1_close = Data_r1_close.nc13;
nc14_r1_close = Data_r1_close.nc14;
Time_r1_close = Data_r1_close.ElapsedTime;
MeanFluo_r1_close = Data_r1_close.MeanVectorAP;
SEFluo_r1_close = Data_r1_close.SEVectorAP;

MeanFluo_individual_r1_close = Data_r1_close.MeanVectorAP_individual;

%% [111] : r3 (female)
nc13_r3 = Data_r3.nc13;
nc14_r3 = Data_r3.nc14;
Time_r3 = Data_r3.ElapsedTime;
MeanFluo_r3 = Data_r3.MeanVectorAP;
SEFluo_r3 = Data_r3.SEVectorAP;

MeanFluo_individual_r3 = Data_r3.MeanVectorAP_individual;
%% Plots for time traces in nc13 (nc14) of [100]
tWindow = nc13_r1:nc14_r1;
%mkdir([DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'])

for APbin = 8:25 % (APbin-1)*2.5 % of Embryo length (15%-60%)
    clf
    hold on

    errorbar(Time_r1(tWindow)-Time_r1(nc13_r1),...
             MeanFluo_r1(tWindow,APbin), SEFluo_r1(tWindow,APbin),...
             'LineWidth',2,'Color',ColorChoice(2,:))
    
    xlim([0 18])
    ylim([0 1000])
    xlabel('time into NC13(min)')
    ylabel('mean spot fluorescence (AU)')
    legend('[100]')
    StandardFigure(gcf,gca)
    %pause
    
    figPath = [DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'];
    saveas(gcf,[figPath,filesep,'averaged_MS2_TimeTraces_[100]' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.tif']); 
    saveas(gcf,[figPath,filesep,'averaged_MS2_TimeTraces_[100]' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.pdf']);
end

%% Plots for time traces in nc13 (nc14) of [010] : r1-mid
tWindow_r1_mid = nc13_r1_mid:nc14_r1_mid;
%mkdir([DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'])

for APbin = 8:25 % (APbin-1)*2.5 % of Embryo length (15%-60%)
    %clf
    hold on

    errorbar(Time_r1_mid(tWindow_r1_mid)-Time_r1_mid(nc13_r1_mid),...
             MeanFluo_r1_mid(tWindow_r1_mid,APbin), SEFluo_r1_mid(tWindow_r1_mid,APbin),...
             'LineWidth',2,'Color',ColorChoice(6,:))
    
    xlim([0 18])
    ylim([0 1000])
    xlabel('time into NC13(min)')
    ylabel('mean spot fluorescence (AU)')
    legend('[010]')
    StandardFigure(gcf,gca)
    pause
    
    figPath = [DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'];
    %saveas(gcf,[figPath,filesep,'averaged_MS2_TimeTraces_[010]' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.tif']); 
    %saveas(gcf,[figPath,filesep,'averaged_MS2_TimeTraces_[010]' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.pdf']);
end
%% Divide the mean spot fluo at each AP bin with the mean spot fluo at 17.5%
% [000]
denominatorFluo_r0 = MeanFluo_r0(:,8);
% 1 Runt site
denominatorFluo_r1 = MeanFluo_r1(:,8);
denominatorFluo_r1_mid = MeanFluo_r1_mid(:,9); % no data point at 17.5%...
denominatorFluo_r1_close = MeanFluo_r1_close(:,8);
% [111]
denominatorFluo_r3 = MeanFluo_r3(:,8);


FluoRatio_r0 = MeanFluo_r1./denominatorFluo_r1;
FluoRatio_SE_r0 = SEFluo_r1./denominatorFluo_r1;

FluoRatio_r1 = MeanFluo_r1./denominatorFluo_r1;
FluoRatio_SE_r1 = SEFluo_r1./denominatorFluo_r1;

FluoRatio_r1_mid = MeanFluo_r1_mid./denominatorFluo_r1_mid;
FluoRatio_SE_r1_mid = SEFluo_r1_mid./denominatorFluo_r1_mid;

FluoRatio_r1_close = MeanFluo_r1_close./denominatorFluo_r1_close;
FluoRatio_SE_r1_close = SEFluo_r1_close./denominatorFluo_r1_close;

FluoRatio_r3 = MeanFluo_r3./denominatorFluo_r3;
FluoRatio_SE_r3 = SEFluo_r3./denominatorFluo_r3;
%% Plots for the ratio(Mean Fluo(:,APbin) / Mean Fluo(:,17.5%) in nc13 for [100]
tBegin = 5 % ~2mins 
tWindow_r1 = nc13_r1+tBegin:nc14_r1;
tWindow_r1_mid = nc13_r1_mid+tBegin:nc14_r1_mid;
tWindow_r1_close = nc13_r1_close+tBegin:nc14_r1_close;
%mkdir([DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'])

for APbin = 8:25 % (APbin-1)*2.5 % of Embryo length (15%-60%)
    clf
    hold on

    errorbar(Time_r1(tWindow_r1)-Time_r1(nc13_r1),...
             FluoRatio_r1(tWindow_r1,APbin), FluoRatio_SE_r1(tWindow_r1,APbin),...
             'LineWidth',2,'Color',ColorChoice(2,:))
         
    errorbar(Time_r1_mid(tWindow_r1_mid)-Time_r1_mid(nc13_r1_mid),...
             FluoRatio_r1_mid(tWindow_r1_mid,APbin), FluoRatio_SE_r1_mid(tWindow_r1_mid,APbin),...
             'LineWidth',2,'Color',ColorChoice(6,:))

    errorbar(Time_r1_close(tWindow_r1_close)-Time_r1_close(nc13_r1_close),...
             FluoRatio_r1_close(tWindow_r1_close,APbin), FluoRatio_SE_r1_close(tWindow_r1_close,APbin),...
             'LineWidth',2,'Color',ColorChoice(5,:))
         
    xlim([0 18])
    ylim([0 2.5])
    xlabel('time into NC13(min)')
    ylabel('ratio (mean spot fluorescence)')
    legend('[100]','[010]','[001]')
    StandardFigure(gcf,gca)
    %pause
    
    figPath = [DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'];
    %saveas(gcf,[figPath,filesep,'ratio_dividedBy_17p5%_TimeTraces_1RuntSite_diffPositions' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.tif']); 
    %saveas(gcf,[figPath,filesep,'ratio_dividedBy_17p5%_TimeTraces_1RuntSite_diffPositions' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.pdf']); 
end

%% [111] divided with its own @ 17.5%
tBegin = 5 % ~2mins 
tWindow_r3 = nc13_r3+tBegin:nc14_r3;

%mkdir([DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'])

for APbin = 8:25 % (APbin-1)*2.5 % of Embryo length (15%-60%)
    clf
    hold on

    errorbar(Time_r3(tWindow_r3)-Time_r3(nc13_r3),...
             FluoRatio_r3(tWindow_r3,APbin), FluoRatio_SE_r3(tWindow_r3,APbin),...
             'LineWidth',2,'Color',ColorChoice(2,:))
         
    xlim([0 18])
    ylim([0 2.5])
    xlabel('time into NC13(min)')
    ylabel('ratio (mean spot fluorescence)')
    legend('[100]','[010]','[001]')
    StandardFigure(gcf,gca)
    %pause
    
    figPath = [DropboxPath,filesep,'Garcia Lab/Figures/Opposing Gradients/Data/AveragedMS2_TimeTraces_sexed'];
    %saveas(gcf,[figPath,filesep,'ratio_dividedBy_17p5%_TimeTraces_1RuntSite_diffPositions' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.tif']); 
    %saveas(gcf,[figPath,filesep,'ratio_dividedBy_17p5%_TimeTraces_1RuntSite_diffPositions' ,'_AP=',num2str((APbin-1)*2.5),'%_NC13' , '.pdf']); 
end

end