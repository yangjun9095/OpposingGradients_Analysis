% comparison between "dynamic" vs "static" inputs
% May, 2021, Yang Joon Kim
function Supple_compare_dynamic_static_inputs
%% Description 

%% Load the input TF datasets
% load the Bicoid datasets
% nc14, 6 embryos
% BcdData.mat

% load the Runt datasets
% MeanFluo (note which .mat file I saved these values)

% interpolate the Runt data as the frame rate is 1 min / frame for Runt,
% which is 2x slower than Bicoid measurements.

for embryo = 1:7
    for APbin = 1:41
        MeanFluo_interp1(:, APbin, embryo) = interp1(Time, MeanFluo(:, APbin, embryo), 0:0.5:77);        
    end
end
%% variables

%% parameters
Kd = 30; % AU
Kr = 100; % AU
w_bp = 100;
p = 0.001;
R = 300;
w_rp = 0.1;

% define the params as needed by the model function
% params = [Kb, Kr, w_bp, w_rp, p, R];

%% Case1. 6A0R case (considering only Bicoid)
%% generate predictions using a model (dynamic input)

% for the Bicoid profiles from individual embryos, we will generate the
% predicted rate of transcription using the stat mech model, the output
% would be rate(time).

% redefine the params
params = [Kb, w_bp, p, R];

% set the frame length consistent for the shortest frame
frameLength = 56; % need to change this to be automatically updated

for i=1:length(BcdData)
    Bcd_temp = BcdData(i).BcdNC14;
    Bcd_temp = Bcd_temp(1:frameLength,:);
    Bcd_temp = Bcd_temp'; % transpose such that the dimension is AP bins x frames (time)
    for j=1:frameLength % for each frame
        rate_dynamicTF(j,:,i) = model_6A0R_HillModel_V3(params, Bcd_temp(:,j));
    end
end

% calculate the mean and sem of the predicted rate "over embryos"
rate_dynamicTF_mean = nanmean(rate_dynamicTF,3);
rate_dynamicTF_sem = nanstd(rate_dynamicTF,[],3)/sqrt(length(BcdData));

%% (dynamic) calculate the time-averaged rate of transcription initiation over 5-10
% minutes (11th to 21st frames)
frameRange = 11:21;

rate_dynamicTF_mean_tAveraged = nanmean(rate_dynamicTF_mean(frameRange,:));
rate_dynamicTF_sem_tAveraged = nanstd(rate_dynamicTF_mean(frameRange,:))./sqrt(length(BcdData));

%% get the averaged profile
% Use the "Bcd_timeAveraged_5_10min_nc14" field

%% generate predictions using a model (static input)
rate_staticTF = model_6A0R_HillModel_V3(params, Bcd_timeAveraged_5_10min_nc14');

rate_staticTF_SD = model_6A0R_HillModel_V3(params, Bcd_timeAveraged_5_10min_nc14_SD');
rate_staticTF_SEM = rate_staticTF_SD./sqrt(length(BcdData));

%% generate plots of dynamic rate vs static rate

% Pick one AP bin
APbin = 9; % 20%
hold on
% static
errorbar(BcdData(6).TimeNC14, ones(frameLength,1)*rate_staticTF(APbin),...
            ones(frameLength,1)*rate_staticTF_SD(APbin),'color',ColorChoice(4,:), 'linewidth',1)
        
% dynamic
errorbar(BcdData(6).TimeNC14, rate_dynamicTF_mean(:,APbin), ...
            rate_dynamicTF_sem(:,APbin),'color',ColorChoice(1,:), 'linewidth',1)
        
% dynamic (time-averaged)
errorbar(BcdData(6).TimeNC14(11:21), ones(11,1)*rate_dynamicTF_mean_tAveraged(APbin),...
            ones(11,1)*rate_dynamicTF_sem_tAveraged(APbin),'color',ColorChoice(2,:), 'linewidth',1)

xlabel(' time into nc14 (min)')
ylabel(' rate of transcription')

legend('static', 'dynamic', 'dynamic: time-averaged','Location', 'SouthEast')
box on
StandardFigure(gcf,gca)

% save the figure
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\dynamic_vs_static_inputTFs';
saveas(gcf,[FigPath,filesep,'predicted_rate_dynamic_static_over_time', '.pdf']);
%% generate a plot comparing the AP profile (static vs dynamic: time-averaged)

hold on
errorbar(0:0.025:1, rate_staticTF, rate_staticTF_SEM,...
            'color',ColorChoice(4,:), 'linewidth',1)
errorbar(0:0.025:1, rate_dynamicTF_mean_tAveraged, rate_dynamicTF_sem_tAveraged,...
            'color',ColorChoice(2,:), 'linewidth',1)

xlim([0.2 0.8])

xlabel(' time into nc14 (min)')
ylabel(' rate of transcription')

legend('static', 'dynamic: time-averaged','Location', 'NorthEast')
box on
StandardFigure(gcf,gca)

% save the figure
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\dynamic_vs_static_inputTFs';
saveas(gcf,[FigPath,filesep,'predicted_rate_dynamic_static_over_AP', '.pdf']);
%% Case2. 6A1R case (considering both Bicoid and Runt)
% model_6A1R_HillModel_V3_direct

%% generate predictions using a model (dynamic input)

% for the Bicoid profiles from individual embryos, we will generate the
% predicted rate of transcription using the stat mech model, the output
% would be rate(time).

% redefine the params
params = [Kb, Kr, w_bp, w_rp, p, R];

% set the frame length consistent for the shortest frame
frameLength = 56; % need to change this to be automatically updated

embryo = 1; % indexing for the combined embryo for Bicoid and Runt measurements

for i=1:length(BcdData) % Bicoid embryos
    Bcd_temp = BcdData(i).BcdNC14;
    Bcd_temp = Bcd_temp(1:frameLength,:);
    Bcd_temp = Bcd_temp'; % transpose such that the dimension is AP bins x frames (time)
    
    for k = 1:7 % Runt embryos
        Runt_temp = MeanFluo_interp1(1:frameLength,:,k);
        Runt_temp = Runt_temp'; % transpose such that the dimension is AP bins x frames (time)
    
        for j=1:frameLength % for each frame
            rate_dynamicTF(j,:,embryo) = model_6A1R_HillModel_V3_direct(params, [Bcd_temp(:,j), Runt_temp(:,j)]);
        end
        % count the embryo index
        embryo = embryo+1;
    end
    
end

% calculate the mean and sem of the predicted rate "over embryos"
rate_dynamicTF_mean = nanmean(rate_dynamicTF,3);
rate_dynamicTF_SD = nanstd(rate_dynamicTF,[],3);
rate_dynamicTF_sem = nanstd(rate_dynamicTF,[],3)/sqrt(42);

%% (dynamic) calculate the time-averaged rate of transcription initiation over 5-10

% minutes (11th to 21st frames)
frameRange = 11:21;

rate_dynamicTF_mean_tAveraged = nanmean(rate_dynamicTF_mean(frameRange,:));
rate_dynamicTF_SD_tAveraged = nanstd(rate_dynamicTF_mean(frameRange,:));
rate_dynamicTF_sem_tAveraged = nanstd(rate_dynamicTF_mean(frameRange,:))./sqrt(length(BcdData));

%% get the averaged profile
% Use the "Bcd_timeAveraged_5_10min_nc14" field

%% generate predictions using a model (static input)
rate_staticTF = model_6A1R_HillModel_V3_direct(params, [Bcd_timeAveraged_5_10min_nc14', RuntFluo']);

rate_staticTF_SD = model_6A1R_HillModel_V3_direct(params, [Bcd_timeAveraged_5_10min_nc14_SD', SDFluo_tAveraged_mixed]);
rate_staticTF_SEM = rate_staticTF_SD./sqrt(42);

%% generate plots of dynamic rate vs static rate

% Pick one AP bin
APbin = 13; % 20%
hold on
% static
errorbar(BcdData(6).TimeNC14, ones(frameLength,1)*rate_staticTF(APbin),...
            ones(frameLength,1)*rate_staticTF_SD(APbin),'color',ColorChoice(4,:), 'linewidth',1)
        
% dynamic
errorbar(BcdData(6).TimeNC14, rate_dynamicTF_mean(:,APbin), ...
            rate_dynamicTF_sem(:,APbin),'color',ColorChoice(1,:), 'linewidth',1)
        
% dynamic (time-averaged)
errorbar(BcdData(6).TimeNC14(11:21), ones(11,1)*rate_dynamicTF_mean_tAveraged(APbin),...
            ones(11,1)*rate_dynamicTF_SD_tAveraged(APbin),'color',ColorChoice(2,:), 'linewidth',1)

ylim([0 12])
yticks([0 2 4 6 8 10 12])
        
xlabel(' time into nc14 (min)')
ylabel(' rate of transcription')

legend('static', 'dynamic', 'dynamic: time-averaged','Location', 'SouthEast')
box on
StandardFigure(gcf,gca)

% save the figure
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\dynamic_vs_static_inputTFs';
saveas(gcf,[FigPath,filesep,'6A1R_predicted_rate_dynamic_static_over_time_30%', '.pdf']);

%% generate a plot comparing the AP profile (static vs dynamic: time-averaged)

hold on
errorbar(0:0.025:1, rate_staticTF, rate_staticTF_SD,...
            'color',ColorChoice(4,:), 'linewidth',1)
errorbar(0:0.025:1, rate_dynamicTF_mean_tAveraged, rate_dynamicTF_SD_tAveraged,...
            'color',ColorChoice(2,:), 'linewidth',1)

xlim([0.2 0.8])
ylim([0 15])

xlabel(' time into nc14 (min)')
ylabel(' rate of transcription')

legend('static', 'dynamic: time-averaged','Location', 'NorthEast')
box on
StandardFigure(gcf,gca)

% save the figure
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\dynamic_vs_static_inputTFs';
saveas(gcf,[FigPath,filesep,'6A1R_predicted_rate_dynamic_static_over_AP', '.pdf']);
end