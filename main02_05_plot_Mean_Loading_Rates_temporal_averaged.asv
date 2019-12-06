function main02_05_plot_Mean_Loading_Rates_temporal_averaged
%% DESCRIPTION
% Written by Yang Joon Kim (yjkim90@berkeley.edu)
% Last edited : 7/31/2019

% The goal is to calculate the mean loading rates (averaged over some time)
% over the AP axis. (This temporal-averaged loading rate is proportional to
% the accumulated mRNA from an ON nucleus, thus need to multiply the
% Fraction ON to recap the total, accumulated mRNA of an embryo.

% Caveats : I'll use the processed (synchronized to the beginning of the
% NC13) datasets of r0, 1, 2, 3. I can always change the datasets that I
% use, for example, the new r0, instead of the old r0 datasets.

%% Load the datasets
filePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

r0Data = load([filePath,filesep,'r0.mat']);
r1Data = load([filePath,filesep,'r1-new-female.mat']);
r2Data = load([filePath,filesep,'r2-new-female.mat']);
r3Data = load([filePath,filesep,'r3-new-female.mat']);

%% Calculate the time-averaged (mean) spot fluo
% 1) NC13
time_averaged_loading_rate_r0_NC13 = nanmean(r0Data.MeanVectorAP_individual(r0Data.nc13:r0Data.nc14,:,:),1);
time_averaged_loading_rate_r1_NC13 = nanmean(r1Data.MeanVectorAP_individual(r1Data.nc13:r1Data.nc14,:,:),1);
time_averaged_loading_rate_r2_NC13 = nanmean(r2Data.MeanVectorAP_individual(r2Data.nc13:r2Data.nc14,:,:),1);
time_averaged_loading_rate_r3_NC13 = nanmean(r3Data.MeanVectorAP_individual(r3Data.nc13:r3Data.nc14,:,:),1);

% 2) NC14
% For now, I'll just take into account of the beginning 30 minutes, since
% all of the constructs seem to exist until that point. Plus, the mRNA made
% after that point would be smaller (negligible?) than the previous time
% points.
Tmax = 30; % min
% Another assumption is that the time interval is same for all datasets
% (which is true).
time_interval = median(diff(r0Data.ElapsedTime));
T_max_index = ceil(Tmax/time_interval); % index for the 30 minutes time vector.

time_averaged_loading_rate_r0_NC14 = nanmean(r0Data.MeanVectorAP_individual(r0Data.nc14:r0Data.nc14 + T_max_index,:,:),1);
time_averaged_loading_rate_r1_NC14 = nanmean(r1Data.MeanVectorAP_individual(r1Data.nc14:r1Data.nc14 + T_max_index,:,:),1);
time_averaged_loading_rate_r2_NC14 = nanmean(r2Data.MeanVectorAP_individual(r2Data.nc14:r2Data.nc14 + T_max_index,:,:),1);
time_averaged_loading_rate_r3_NC14 = nanmean(r3Data.MeanVectorAP_individual(r3Data.nc14:r3Data.nc14 + T_max_index,:,:),1);
%% Quick plot for a sanity check
% number of embryos
numEmbryos_r0 = length(time_averaged_loading_rate_r0_NC13(1,1,:));
numEmbryos_r1 = length(time_averaged_loading_rate_r1_NC13(1,1,:));
numEmbryos_r2 = length(time_averaged_loading_rate_r2_NC13(1,1,:));
numEmbryos_r3 = length(time_averaged_loading_rate_r3_NC13(1,1,:));

hold on
for i=1:numEmbryos_r0
    plot(0:0.025:1, time_averaged_loading_rate_r0_NC13(1,:,i))
end

hold on
for i=1:numEmbryos_r1
    plot(0:0.025:1, time_averaged_loading_rate_r1_NC13(1,:,i))
end

hold on
for i=1:numEmbryos_r2
    plot(0:0.025:1, time_averaged_loading_rate_r2_NC13(1,:,i))
end

hold on
for i=1:numEmbryos_r3
    plot(0:0.025:1, time_averaged_loading_rate_r3_NC13(1,:,i))
end


%% Average over multiple embryos
%% First, the zeroes should be converted to NaNs. This should actually be
% done in the AverageDatasets. Somehow the temporal averaging is making
% these false zeroes...
% 1) NC13
time_averaged_loading_rate_r0_NC13(time_averaged_loading_rate_r0_NC13==0) = nan;
time_averaged_loading_rate_r1_NC13(time_averaged_loading_rate_r1_NC13==0) = nan;
time_averaged_loading_rate_r2_NC13(time_averaged_loading_rate_r2_NC13==0) = nan;
time_averaged_loading_rate_r3_NC13(time_averaged_loading_rate_r3_NC13==0) = nan;
% 2) NC14
time_averaged_loading_rate_r0_NC14(time_averaged_loading_rate_r0_NC14==0) = nan;
time_averaged_loading_rate_r1_NC14(time_averaged_loading_rate_r1_NC14==0) = nan;
time_averaged_loading_rate_r2_NC14(time_averaged_loading_rate_r2_NC14==0) = nan;
time_averaged_loading_rate_r3_NC14(time_averaged_loading_rate_r3_NC14==0) = nan;

%% Second, use nanmean and nanstd to get the average over multiple embryos.
% Actually, this should be more sophisticated, like having 
% r0, NC13
Mean_time_averaged_loading_rate_r0_NC13 = nanmean(time_averaged_loading_rate_r0_NC13,3);
SEM_time_averaged_loading_rate_r0_NC13 = nanstd(time_averaged_loading_rate_r0_NC13,[],3)./sqrt(numEmbryos_r0);
% r0, NC14
Mean_time_averaged_loading_rate_r0_NC14 = nanmean(time_averaged_loading_rate_r0_NC14,3);
SEM_time_averaged_loading_rate_r0_NC14 = nanstd(time_averaged_loading_rate_r0_NC14,[],3)./sqrt(numEmbryos_r0);

% r1
Mean_time_averaged_loading_rate_r1_NC13 = nanmean(time_averaged_loading_rate_r1_NC13,3);
SEM_time_averaged_loading_rate_r1_NC13 = nanstd(time_averaged_loading_rate_r1_NC13,[],3)./sqrt(numEmbryos_r1);
% r1, NC14
Mean_time_averaged_loading_rate_r1_NC14 = nanmean(time_averaged_loading_rate_r1_NC14,3);
SEM_time_averaged_loading_rate_r1_NC14 = nanstd(time_averaged_loading_rate_r1_NC14,[],3)./sqrt(numEmbryos_r1);

% r2, NC13
Mean_time_averaged_loading_rate_r2_NC13 = nanmean(time_averaged_loading_rate_r2_NC13,3);
SEM_time_averaged_loading_rate_r2_NC13 = nanstd(time_averaged_loading_rate_r2_NC13,[],3)./sqrt(numEmbryos_r2);
% r2, NC14
Mean_time_averaged_loading_rate_r2_NC14 = nanmean(time_averaged_loading_rate_r2_NC14,3);
SEM_time_averaged_loading_rate_r2_NC14 = nanstd(time_averaged_loading_rate_r2_NC14,[],3)./sqrt(numEmbryos_r2);

% r3, NC13
Mean_time_averaged_loading_rate_r3_NC13 = nanmean(time_averaged_loading_rate_r3_NC13,3);
SEM_time_averaged_loading_rate_r3_NC13 = nanstd(time_averaged_loading_rate_r3_NC13,[],3)./sqrt(numEmbryos_r3);
% r3, NC14
Mean_time_averaged_loading_rate_r3_NC14 = nanmean(time_averaged_loading_rate_r3_NC14,3);
SEM_time_averaged_loading_rate_r3_NC14 = nanstd(time_averaged_loading_rate_r3_NC14,[],3)./sqrt(numEmbryos_r3);
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
%% Plot time-averaged loading rates (also averaged over multiple embryos)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Time_averaged_LoadingRate';
APaxis = 0:0.025:1;

% NC13
Fig_Time_averaged_LoadingRate_NC13 = figure;
hold on
errorbar(APaxis, Mean_time_averaged_loading_rate_r0_NC13, SEM_time_averaged_loading_rate_r0_NC13,'Color',ColorChoice(1,:))
errorbar(APaxis, Mean_time_averaged_loading_rate_r1_NC13, SEM_time_averaged_loading_rate_r1_NC13,'Color',ColorChoice(2,:))
errorbar(APaxis, Mean_time_averaged_loading_rate_r2_NC13, SEM_time_averaged_loading_rate_r2_NC13,'Color',ColorChoice(3,:))
errorbar(APaxis, Mean_time_averaged_loading_rate_r3_NC13, SEM_time_averaged_loading_rate_r3_NC13,'Color',ColorChoice(4,:))

legend('r0','r1','r2','r3')%,'r1-male','r2-male','r3-male') 
xlabel('AP Position')
ylabel('RNAP loading rate (AU/min)')

title('time-averaged RNAP loading rates along AP axis, during NC 13')
StandardFigure(Fig_Time_averaged_LoadingRate_NC13,Fig_Time_averaged_LoadingRate_NC13.CurrentAxes)

% Save figures
saveas(Fig_Time_averaged_LoadingRate_NC13,[FigPath,filesep, 'Time_Averaged_RNAP_loading_rates_r0123' , '_NC13' , '.tif']); 
saveas(Fig_Time_averaged_LoadingRate_NC13,[FigPath,filesep, 'Time_Averaged_RNAP_loading_rates_r0123' , '_NC13' , '.pdf']); 

% NC14
Fig_Time_averaged_LoadingRate_NC14 = figure;
hold on
errorbar(APaxis, Mean_time_averaged_loading_rate_r0_NC14, SEM_time_averaged_loading_rate_r0_NC14,'Color',ColorChoice(1,:))
errorbar(APaxis, Mean_time_averaged_loading_rate_r1_NC14, SEM_time_averaged_loading_rate_r1_NC14,'Color',ColorChoice(2,:))
errorbar(APaxis, Mean_time_averaged_loading_rate_r2_NC14, SEM_time_averaged_loading_rate_r2_NC14,'Color',ColorChoice(3,:))
errorbar(APaxis, Mean_time_averaged_loading_rate_r3_NC14, SEM_time_averaged_loading_rate_r3_NC14,'Color',ColorChoice(4,:))

legend('r0','r1','r2','r3')%,'r1-male','r2-male','r3-male') 
xlabel('AP Position')
ylabel('RNAP loading rate (AU/min)')

title('time-averaged RNAP loading rates along AP axis, during NC 14')
StandardFigure(Fig_Time_averaged_LoadingRate_NC14,Fig_Time_averaged_LoadingRate_NC14.CurrentAxes)

% Save figures
saveas(Fig_Time_averaged_LoadingRate_NC14,[FigPath,filesep, 'Time_Averaged_RNAP_loading_rates_r0123' , '_NC14' , '.tif']); 
saveas(Fig_Time_averaged_LoadingRate_NC14,[FigPath,filesep, 'Time_Averaged_RNAP_loading_rates_r0123' , '_NC14' , '.pdf']); 
end