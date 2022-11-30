function main08_04_Compare_averaged_Txn_output_sexes(DataType,varargin)
% Compare gene expression profile (transcription) for the same constructs
% for different sexes

%DataType = 'r0-new'

DataType_male = [DataType,'-male'];
DataType_female = [DataType,'-female'];

filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\TxnOutput_sexed';

Data_male = load([filePath,filesep,DataType_male,'.mat']);
Data_female = load([filePath,filesep,DataType_female,'.mat']);


%% Compare the Mean spot fluorescence over Time

MeanSpotFluo_male = Data_male.MeanVectorAP;
SESpotFluo_male = Data_male.SEVectorAP;
NC13_male = Data_male.nc13;
NC14_male = Data_male.nc14;
Time_male = Data_male.ElapsedTime;

MeanSpotFluo_female = Data_female.MeanVectorAP;
SESpotFluo_female = Data_female.SEVectorAP;
NC13_female = Data_female.nc13;
NC14_female = Data_female.nc14;
Time_female = Data_female.ElapsedTime;


%% Post-processing
% This is weird, but especially r3, some APbins have time frames with
% negative MeanVectorAP... I need to go back, and check the segmentation,
% etc. Also, note that this spot fluo is calculated with flat plane offset.
% (6/6/2019). I can try with AR's new offset calculation.
% For now, I'll assign NaNs for negative MeanVectorAP.


%% Plot

% Define the NC to plot
NC=14;
if NC==13
    tRange_male = NC13_male:NC14_male;
    tRange_female = NC13_female:NC14_female;
    t0_male = 0;
    t0_female = 0;
elseif NC==14
    tRange_male = NC14_male:length(Time_male);
    tRange_female = NC14_female:length(Time_female);
    t0_male = Time_male(NC14_male);
    t0_female = Time_female(NC14_female);
end

% For all APbins, plot the Mean spot fluo over that NC.
APstart = 10; %22.5%
APend = 21; % 50%

for APbin = APstart:APend %9:17 % 20-40% of AP axis
    AP = (APbin-1)*2.5;
    % Check if both male and female data exist (as not NaNs) in that APbin
    if sum(~isnan(MeanSpotFluo_male(:,APbin))) > 0 &&...
            sum(~isnan(MeanSpotFluo_female(:,APbin))) > 0 

        Particle_trace_figure(APbin-APstart+1) = figure(APbin-APstart+1)
        hold on
        errorbar(Time_male(tRange_male) - t0_male, MeanSpotFluo_male(tRange_male,APbin), ...
                    SESpotFluo_male(tRange_male,APbin))

        errorbar(Time_female(tRange_female) - t0_female, MeanSpotFluo_female(tRange_female,APbin), ...
                    SESpotFluo_female(tRange_female,APbin))

        title({'Mean MS2 Spot fluorescence over time';[' AP = ',num2str(AP),'%']})
        xlabel('Time (min)')
        ylabel('Mean spot fluorescence (AU)')
        legend('male','female')
        %legend([h(k) H.mainLine],'single','Mean')
        %l = l+1;
        xlim([0 max(tRange_male)])

    else
        display(['APbin ',num2str(AP),'% does not have non-Nan values'])
    end
end
       
hold off
%% Make the figure look better & save

% Generate folder to save the data
DataFolder=['E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\hbP2-r0123-AveragedMS2Traces_sex\',DataType(1:end-1),'\Mean_spot_fluo_NC',num2str(NC),'_SEM'];
mkdir(DataFolder)
FigPath = DataFolder;

for i=1:length(Particle_trace_figure)
    StandardFigure(Particle_trace_figure(i),Particle_trace_figure(i).CurrentAxes)

    % saveas(Particle_trace_figure(i),[FigPath,filesep, 'Embryos_together','_single_vs_averaged_MS2_traces at ' ,num2str((i+8-1)*2.5), '%_NC',num2str(NC) , '.tif']); 
    % saveas(Particle_trace_figure(i),[FigPath,filesep, 'Embryos_together','_single_vs_averaged_MS2_traces at ' ,num2str((i+8-1)*2.5), '%_NC',num2str(NC) , '.pdf']); 
    saveas(Particle_trace_figure(i),[FigPath,filesep,'Mean_MS2_traces at ' ,num2str((i+APstart-2)*2.5), '%_NC',num2str(NC) , '.tif']); 
    saveas(Particle_trace_figure(i),[FigPath,filesep,'Mean_MS2_traces at ' ,num2str((i+APstart-2)*2.5), '%_NC',num2str(NC) , '.pdf']); 

end

pause
close all
%% Accumulated mRNA
AccumulatedmRNA_male = Data_male.AccumulatedmRNA;
AccumulatedmRNA_FractionON_male = Data_male.AccumulatedmRNA_FractionON;
AccumulatedmRNA_FractionON_SD_male = Data_male.AccumulatedmRNA_FractionON_SD;

AccumulatedmRNA_female = Data_female.AccumulatedmRNA;
AccumulatedmRNA_FractionON_female = Data_female.AccumulatedmRNA_FractionON;
AccumulatedmRNA_FractionON_SD_female = Data_female.AccumulatedmRNA_FractionON_SD;

%% Plot - Accumulated mRNA

% Make directory to save the results
DataFolder=['E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Transcription-Output\hbP2-r0123-AccumulatedmRNA_sex'];
mkdir(DataFolder)
FigPath = DataFolder;

% Define the AP axis
APaxis = 0:0.025:1;

% 1) end of NC13 (beginning of NC14)
AccumulatedmRNA_nc13_figure = figure
hold on
errorbar(APaxis, AccumulatedmRNA_FractionON_male(NC14_male,:),...
                    AccumulatedmRNA_FractionON_SD_male(NC14_male,:))
                
errorbar(APaxis, AccumulatedmRNA_FractionON_female(NC14_female,:),...
                    AccumulatedmRNA_FractionON_SD_female(NC14_female,:))
                
title('Accumulated mRNA _ Fraction ON')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('Male','Female')

StandardFigure(AccumulatedmRNA_nc13_figure,AccumulatedmRNA_nc13_figure. CurrentAxes)
saveas(AccumulatedmRNA_nc13_figure,[FigPath,filesep,'AccumulatedmRNA_',DataType(1:end-1),'_NC13','_SD','.tif'])
saveas(AccumulatedmRNA_nc13_figure,[FigPath,filesep,'AccumulatedmRNA_',DataType(1:end-1),'_NC13','_SD','.pdf'])

% 2) At the end of NC14
AccumulatedmRNA_nc14_figure = figure
hold on
errorbar(APaxis, AccumulatedmRNA_FractionON_male(end,:),...
                    AccumulatedmRNA_FractionON_SD_male(end,:))
                
errorbar(APaxis, AccumulatedmRNA_FractionON_female(end,:),...
                    AccumulatedmRNA_FractionON_SD_female(end,:))
                
title('Accumulated mRNA _ Fraction ON')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('Male','Female')

StandardFigure(AccumulatedmRNA_nc14_figure,AccumulatedmRNA_nc14_figure. CurrentAxes)
saveas(AccumulatedmRNA_nc14_figure,[FigPath,filesep,'AccumulatedmRNA_',DataType(1:end-1),'_NC14','_SD','.tif'])
saveas(AccumulatedmRNA_nc14_figure,[FigPath,filesep,'AccumulatedmRNA_',DataType(1:end-1),'_NC14','_SD','.pdf'])

pause
close all
end