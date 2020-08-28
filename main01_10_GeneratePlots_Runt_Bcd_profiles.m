function main01_10_GeneratePlots_Runt_Bcd_profiles
% this script is for generating plots of processed, time-averaged (or not)
% spatio-temporal profiles of Runt and Bcd. Those datasets are pretty much
% directly processed from main01_01,02,03,04, etc.

%% Part1. Deprecated
%% Runt profile (NC14)
% 
% % load datasets (processed, and time-averaged)
% 
% % Bicoid : Time-averaging is needed, but we know the Bcd gradient shape
% % (length constant) is conserved, and only the amplitude changes over time
% % (from Paul's fitting with Exponential), thus we'll just use the Bcd
% % profile at one time point. The K_a will scale things properly anyway.
% FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
% 
% % Now, I will start from the NC14.
% % Bcd and Runt should be time-averaged for 0-10 minutes in NC14 (the initial
% % rise phase) : YJK on 8/8/2019 : This assumption could change...
% 
% % 1) Bcd :Take the Bicoid data from Liz&Jonathan
% % time resolution : 30 sec.
% BcdAnt = load([FilePath, filesep, 'BcdGFPAnt.mat']);
% 
% BcdNC14 = BcdAnt.DataBcd.nc14;
% % Take temporal average over 0-10 minutes into NC13.
% BcdFluo = nanmean(BcdAnt.DataBcd.MeanVectorAP(BcdNC14: BcdNC14 + 20,:));
% BcdFluoSD = nanmean(BcdAnt.DataBcd.SDVectorAP(BcdNC14: BcdNC14 + 20,:));
% 
% Runt_female = load([FilePath, filesep, 'Runt_time_averaged_female.mat']);
% Runt_male = load([FilePath, filesep, 'Runt_time_averaged_male.mat']);
% 
% % Runt is already time-averaged. For 0-10 minutes, it's in the 5th cell.
% RuntFluo_female = Runt_female.Averaged_Fluo(5,:);
% RuntFluoSE_female = Runt_female.SE_Fluo(5,:);
% 
% RuntFluo_male = Runt_male.Averaged_Fluo(5,:);
% RuntFluoSE_male = Runt_male.SE_Fluo(5,:);
% 
% % Smoothening
% RuntFluo_female = movmean(RuntFluo_female,3); % smoothening with 3 AP bins
% RuntFluo_male = movmean(RuntFluo_male,3); % smoothening with 3 AP bins
% 
% % Extrapoloate the Runt data in anterior bins (since I need from 20%)
% % female
% RuntFluo_female_extrap = interp1([0.3:0.025:0.625], RuntFluo_female(13:26), [0.2:0.025:0.275], 'pchip', 'extrap')
% RuntFluo_female_extrapolated = nan(41,1);
% RuntFluo_female_extrapolated(9:12) = RuntFluo_female_extrap; % 20%-27.5%
% RuntFluo_female_extrapolated(13:27) = RuntFluo(13:27); % 30-65%
% % male
% RuntFluo_male_extrap = interp1([0.3:0.025:0.625], RuntFluo_male(13:26), [0.2:0.025:0.275], 'pchip', 'extrap')
% RuntFluo_male_extrapolated = nan(41,1);
% RuntFluo_male_extrapolated(9:12) = RuntFluo_male_extrap; % 20%-27.5%
% RuntFluo_male_extrapolated(13:27) = RuntFluo(13:27); % 30-65%
% 
% % plot to check
% APaxis = 0:0.025:1;
% RuntScale = 5; % Scaling for plot
% 
% hold on
% yyaxis left
% errorbar(APaxis, BcdFluo, BcdFluoSD)
% ylabel('Bicoid concentration (AU)')
% 
% yyaxis right
% errorbar(APaxis, RuntFluo_female * RuntScale, RuntFluoSE_female)
% errorbar(APaxis, RuntFluo_male * RuntScale, RuntFluoSE_male)
% 
% plot(APaxis, RuntFluo_female_extrapolated * RuntScale)
% plot(APaxis, RuntFluo_male_extrapolated * RuntScale)
% ylabel('Runt concentration (AU)')
% 
% xlim([0.2 0.6])
% 
% title({'Bcd and Runt profile over AP';'(time-averaged, 0-10 min in NC14)'})
% xlabel('AP (Embryo Length)')
% legend('Bcd','Runt-female','','Runt-male','')
% 
% StandardFigure(gcf,gca)
% 
% % Save figures
% % FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Fitting_InitialSlope_Asymmetric';
% % saveas(gcf,[FigPath,filesep,'InputTF_Time_Averaged_NC14_0-10min_BcdRunt' , '_NC14' , '.tif']); 
% % saveas(gcf,[FigPath,filesep,'InputTF_Time_Averaged_NC14_0-10min_BcdRunt' , '_NC14' , '.pdf']); 
% %% Save the extrapolated Runt profile
% Runt = load([FilePath,filesep,'Runt_time_averaged_female.mat'])
% Runt.Extrapolated_Fluo_0_10min_NC14 = RuntFluo_extrapolated;
% 
% save([FilePath,filesep,'Runt_time_averaged_female.mat'],...
%    '-struct', 'Runt','-v7.3');

%% Part2. generate plots of Bcd and Runt (time-traces)
% For a figure of Bcd and Runt's time trace
filePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

% Now, I will start from the NC14.
% Bcd and Runt should be time-averaged for 0-10 minutes in NC14 (the initial
% rise phase) : YJK on 8/8/2019 : This assumption could change...
% But I don't think it matters a lot, since the shape (spatial profile) of
% the gradient would not change for different time windows, since they
% scale proportionately.

% % 1) Bcd :Take the Bicoid data from Liz&Jonathan
% % time resolution : 30 sec.
% BcdAnt = load([FilePath, filesep, 'BcdGFPAnt.mat']);
% 2) YJK + Paul's datasets
BcdData = load([filePath, filesep,'Bcd-Averaged.mat'])

%% Extract useful fields
BcdTime = BcdData.ElapsedTime;
BcdNC13 = BcdData.nc13;
BcdNC14 = BcdData.nc14;

BcdFluo = BcdData.MeanVectorAP;
BcdFluoSD = BcdData.SDVectorAP;


end