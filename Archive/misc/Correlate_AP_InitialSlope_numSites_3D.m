function Correlate_AP_InitialSlope_numSites_3D
%% DESCRIPTION
% Generate exploratory 3D plots of
% x : AP, y : initial slope (either NC13 or NC14), and z : number of Runt
% binding sites.
%% Load the datasets (initial slope, fitted either by Asymmetric or just manual fitting, etc.)
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
Data = load([filePath, filesep, 'AveragedInitialRate_trapezoidalfit.mat']);
%% Extract useful fields from the Data
initialFit_r0 = Data.average_fittedRate_r0;
initialFit_r1 = Data.average_fittedRate_r1_female;
initialFit_r2 = Data.average_fittedRate_r2_female;
initialFit_r3 = Data.average_fittedRate_r3_female;

initialFit_SEM_r0 = Data.SEM_fittedRate_r0;
initialFit_SEM_r1 = Data.SEM_fittedRate_r1_female;
initialFit_SEM_r2 = Data.SEM_fittedRate_r2_female;
initialFit_SEM_r3 = Data.SEM_fittedRate_r3_female;

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow; colorDict.red; colorDict.brown]; % 4 embryos max. it could be extended easily
%% Generate 3D plots
APaxis = 0:0.025:1;
NC=13;
hold on
for i=1:41 % AP bins
    plot3(APaxis(i),0,initialFit_r0(i,NC-11),'o','Color',ColorChoice(1,:))
    plot3(APaxis(i),0,initialFit_r0(i,NC-11) - initialFit_SEM_r0(i,NC-11):...
            initialFit_r0(i,NC-11) + initialFit_SEM_r0(i,NC-11),'-','LineWidth',2,'Color',ColorChoice(1,:))
    plot3(APaxis(i),1,initialFit_r1(i,NC-11),'o','Color',ColorChoice(2,:))
    plot3(APaxis(i),2,initialFit_r2(i,NC-11),'o','Color',ColorChoice(3,:))
    plot3(APaxis(i),3,initialFit_r3(i,NC-11),'o','Color',ColorChoice(4,:))
end

view(20, 60)
title('AP vs number of Runt sites vs initial rate of Txn')
xlabel('AP (Embryo Length)')
ylabel('Number of Runt sites')
zlabel('Initial rate of Txn (AU)')
%legend('','','','')
StandardFigure(gcf,gca)

%% Color bar
%% Plot # of Runt sites vs Rate (@ all AP bins)
NC=13;
hold on
for i=1:41
%     errorbar(0,initialFit_r0(i,NC-11),initialFit_SEM_r0(i,NC-11),'Color',ColorChoice(1,:))
%     errorbar(1,initialFit_r1(i,NC-11),initialFit_SEM_r1(i,NC-11),'Color',ColorChoice(2,:))
%     errorbar(2,initialFit_r2(i,NC-11),initialFit_SEM_r2(i,NC-11),'Color',ColorChoice(3,:))
%     errorbar(3,initialFit_r3(i,NC-11),initialFit_SEM_r3(i,NC-11),'Color',ColorChoice(4,:))
    errorbar([0, 1, 2, 3], [initialFit_r0(i,NC-11),initialFit_r1(i,NC-11),initialFit_r2(i,NC-11),initialFit_r3(i,NC-11) ],...
                [initialFit_SEM_r0(i,NC-11),initialFit_SEM_r1(i,NC-11),initialFit_SEM_r2(i,NC-11),initialFit_SEM_r3(i,NC-11)])%,'Color',ColorChoice(1:4,:))
end