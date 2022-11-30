%% Initial loading rates _ all constructs _ mixed sex
% This is for generating plots of initial loading rates with embryos of
% mixed sex.
% Last updated : October, 2019. For YJK's group meeting.

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

ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
                colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
                colorDict.darkgreen; colorDict.blue; colorDict.cyan; colorDict.darkgreen]; 

%% Load all the datasets that are fitted (by MeanFitAsymmetric.m)
% The caveats are 
% 1) I'll mix males/females, since Runt seems to be dosage compensated in
% early NC14.

% Note that this r0 is from the old construct. Probably wouldn't make much
% difference, but worth making notes for the final version.
Data_r0 = LoadMS2Sets('r0','dontCompare');

Data_r1 = LoadMS2Sets('r1-new','dontCompare');
Data_r2 = LoadMS2Sets('r2-new','dontCompare');
Data_r3 = LoadMS2Sets('r3-new','dontCompare');

% Variants (different binding sites configurations)
Data_r1_mid = LoadMS2Sets('r1-mid','dontCompare');
Data_r1_close = LoadMS2Sets('r1-close','dontCompare');
Data_r2_close = LoadMS2Sets('r2_1+2','dontCompare');
Data_r2_far = LoadMS2Sets('r2_1+3','dontCompare');

% Variant (1bp mutation in Runt binding sites)
Data_r3prime = LoadMS2Sets('r3prime','dontCompare')

%% Extract the fitted initial loading rates
% Original constructs
[fittedRate_r0,fittedRateSD_r0,fittedTon_r0] = Extract_Fields_MeanFits(Data_r0,'Asymmetric');
[fittedRate_r1,fittedRateSD_r1,fittedTon_r1] = Extract_Fields_MeanFits(Data_r1,'Asymmetric');
[fittedRate_r2,fittedRateSD_r2,fittedTon_r2] = Extract_Fields_MeanFits(Data_r2,'Asymmetric');
[fittedRate_r3,fittedRateSD_r3,fittedTon_r3] = Extract_Fields_MeanFits(Data_r3,'Asymmetric');

% Variants
% 1 binding site
[fittedRate_r1_mid,fittedRateSD_r1_mid,fittedTon_r1_mid] = Extract_Fields_MeanFits(Data_r1_mid,'Asymmetric');
[fittedRate_r1_close,fittedRateSD_r1_close,fittedTon_r1_close] = Extract_Fields_MeanFits(Data_r1_close,'Asymmetric');
% 2 binding sites
[fittedRate_r2_close,fittedRateSD_r2_close,fittedTon_r2_close] = Extract_Fields_MeanFits(Data_r2_close,'Asymmetric');
[fittedRate_r2_far,fittedRateSD_r2_far,fittedTon_r2_far] = Extract_Fields_MeanFits(Data_r2_far,'Asymmetric');

% Variant (Mutated Runt binding sites)
[fittedRate_r3prime,fittedRateSD_r3prime,fittedTon_r3prime] = Extract_Fields_MeanFits(Data_r3prime,'Asymmetric');

%% Averaging over embryos using nanmean, nanstd
% r0
average_fittedRate_r0 = nanmean(fittedRate_r0,3);
SEM_fittedRate_r0 = nanstd(fittedRate_r0,0,3)/sqrt(length(Data_r0));
% r1
average_fittedRate_r1 = nanmean(fittedRate_r1,3);
SEM_fittedRate_r1 = nanstd(fittedRate_r1,0,3)/sqrt(length(Data_r1));
% r2
average_fittedRate_r2 = nanmean(fittedRate_r2,3);
SEM_fittedRate_r2 = nanstd(fittedRate_r2,0,3)/sqrt(length(Data_r2));
% r3
average_fittedRate_r3 = nanmean(fittedRate_r3,3);
SEM_fittedRate_r3 = nanstd(fittedRate_r3,0,3)/sqrt(length(Data_r3));

% r3' (r3prime)
average_fittedRate_r3prime = nanmean(fittedRate_r3prime,3);
SEM_fittedRate_r3prime = nanstd(fittedRate_r3prime,0,3)/sqrt(length(Data_r3prime));

% Variants
% r1-mid
average_fittedRate_r1_mid = nanmean(fittedRate_r1_mid,3);
SEM_fittedRate_r1_mid = nanstd(fittedRate_r1_mid,0,3)/sqrt(length(Data_r1_mid));
% r1-close
average_fittedRate_r1_close = nanmean(fittedRate_r1_close,3);
SEM_fittedRate_r1_close = nanstd(fittedRate_r1_close,0,3)/sqrt(length(Data_r1_close));

% r2-close
average_fittedRate_r2_close = nanmean(fittedRate_r2_close,3);
SEM_fittedRate_r2_close = nanstd(fittedRate_r2_close,0,3)/sqrt(length(Data_r2_close));
% r2-far
average_fittedRate_r2_far = nanmean(fittedRate_r2_far,3);
SEM_fittedRate_r2_far = nanstd(fittedRate_r2_far,0,3)/sqrt(length(Data_r2_far));
%% Plot the initial rate from individual embryos vs mean
APaxis = 0:0.025:1;
NC = 3; % NC14
% Start from r0
figure_r0_ind_mean = figure;
hold on
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC))
for i=1:length(Data_r0)
    errorbar(APaxis, fittedRate_r0(:,NC,i), fittedRateSD_r0(:,NC,i))
end
%ylim([0 400])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','e1','e2','e3','e4')

StandardFigure(figure_r0_ind_mean,figure_r0_ind_mean.CurrentAxes)

% save the plot
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MixedSex_Oct_2019';
saveas(figure_r0_ind_mean,[FigPath,filesep,'InitialSlope_r0_individual_mean' , '.pdf']);
saveas(figure_r0_ind_mean,[FigPath,filesep,'InitialSlope_r0_individual_mean' , '.tif']);
%% r1
figure_r1_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r1(:,NC), SEM_fittedRate_r1(:,NC))
% plot individuals
for i=1:length(Data_r1)
    errorbar(APaxis, fittedRate_r1(:,NC,i), fittedRateSD_r1(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3';'e4';'e5'; 'e6';'e7'];
%ylim([0 400])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','male1','male2','male3',...
        'female1','female2','female3')

StandardFigure(figure_r1_ind_mean,figure_r1_ind_mean.CurrentAxes)

% save the plot
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MxiedSex_Oct_2019';
saveas(figure_r1_ind_mean,[FigPath,filesep,'InitialSlope_r1_individual_mean' , '.pdf']);
saveas(figure_r1_ind_mean,[FigPath,filesep,'InitialSlope_r1_individual_mean' , '.tif']);


%% r3
figure_r3_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC))
% plot individuals
for i=1:length(Data_r3)
    errorbar(APaxis, fittedRate_r3(:,NC,i), fittedRateSD_r3(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3';'e4';'e5'; 'e6';'e7'];
%ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','male1','male2','male3','male4',...
        'female1','female2','female3')

StandardFigure(figure_r3_ind_mean,figure_r3_ind_mean.CurrentAxes)

% save the plot
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MxiedSex_Oct_2019';
saveas(figure_r3_ind_mean,[FigPath,filesep,'InitialSlope_r3_individual_mean' , '.pdf']);
saveas(figure_r3_ind_mean,[FigPath,filesep,'InitialSlope_r3_individual_mean' , '.tif']);

%% r3 prime
figure_r3prime_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r3prime(:,NC), SEM_fittedRate_r3prime(:,NC))
% plot individuals
for i=1:length(Data_r3prime)
    errorbar(APaxis, fittedRate_r3prime(:,NC,i), fittedRateSD_r3prime(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3'];
%ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','embryos1','embryo2','embryo3')

StandardFigure(figure_r3prime_ind_mean,figure_r3prime_ind_mean.CurrentAxes)

% save the plot
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MixedSex_Oct_2019';
saveas(figure_r3prime_ind_mean,[FigPath,filesep,'InitialSlope_r3prime_individual_mean' , '.pdf']);
saveas(figure_r3prime_ind_mean,[FigPath,filesep,'InitialSlope_r3prime_individual_mean' , '.tif']);


%% r1-close
figure_r1_close_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r1_close(:,NC), SEM_fittedRate_r1_close(:,NC))
% plot individuals
for i=1:length(Data_r1_close)
    errorbar(APaxis, fittedRate_r1_close(:,NC,i), fittedRateSD_r1_close(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3';'e4';'e5'; 'e6';'e7'];
%ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','female1','female2','male1','female3')

StandardFigure(figure_r1_close_ind_mean,figure_r1_close_ind_mean.CurrentAxes)

% save the plot
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MxiedSex_Oct_2019';
saveas(figure_r1_close_ind_mean,[FigPath,filesep,'InitialSlope_r1_close_individual_mean' , '.pdf']);
saveas(figure_r1_close_ind_mean,[FigPath,filesep,'InitialSlope_r1_close_individual_mean' , '.tif']);
%% r1-mid
figure_r1_mid_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r1_mid(:,NC), SEM_fittedRate_r1_mid(:,NC))
% plot individuals
for i=1:length(Data_r1_mid)
    errorbar(APaxis, fittedRate_r1_mid(:,NC,i), fittedRateSD_r1_mid(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3';'e4';'e5'; 'e6';'e7'];
%ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','Unknown','male1','female1','Unknown')

StandardFigure(figure_r1_mid_ind_mean,figure_r1_mid_ind_mean.CurrentAxes)

% save the plot
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MxiedSex_Oct_2019';
saveas(figure_r1_mid_ind_mean,[FigPath,filesep,'InitialSlope_r1_mid_individual_mean' , '.pdf']);
saveas(figure_r1_mid_ind_mean,[FigPath,filesep,'InitialSlope_r1_mid_individual_mean' , '.tif']);
%% r2-close (1+2)
figure_r2_close_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r2_close(:,NC), SEM_fittedRate_r2_close(:,NC))
% plot individuals
for i=1:length(Data_r2_close)
    errorbar(APaxis, fittedRate_r2_close(:,NC,i), fittedRateSD_r2_close(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3';'e4';'e5'; 'e6';'e7'];
%ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','Unknown1','Unknown2','Unknown3',...
        'Unknown4')

StandardFigure(figure_r2_close_ind_mean,figure_r2_close_ind_mean.CurrentAxes)

% save the plot
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MxiedSex_Oct_2019';
saveas(figure_r2_close_ind_mean,[FigPath,filesep,'InitialSlope_r2_close_1+2_individual_mean' , '.pdf']);
saveas(figure_r2_close_ind_mean,[FigPath,filesep,'InitialSlope_r2_close_1+2_individual_mean' , '.tif']);

%% r2-far(1+3)
figure_r2_far_ind_mean = figure;
hold on
% plot the Averaged
errorbar(APaxis, average_fittedRate_r2_far(:,NC), SEM_fittedRate_r2_far(:,NC))
% plot individuals
for i=1:length(Data_r2_far)
    errorbar(APaxis, fittedRate_r2_far(:,NC,i), fittedRateSD_r2_far(:,NC,i))
end
LegendLabel = ['e1';'e2';'e3';'e4';'e5'; 'e6';'e7'];
%ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('Average','Unknown1','Unknown2','Unknown3')%,...
        %'Unknown4')

StandardFigure(figure_r2_far_ind_mean,figure_r2_far_ind_mean.CurrentAxes)

% save the plot
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InitialSlope_Asymmetric\MxiedSex_Oct_2019';
saveas(figure_r2_far_ind_mean,[FigPath,filesep,'InitialSlope_r2_far_1+3_individual_mean' , '.pdf']);
saveas(figure_r2_far_ind_mean,[FigPath,filesep,'InitialSlope_r2_far_1+3_individual_mean' , '.tif']);

%% Comparing the 1 binding site constructs
% ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
%                 colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
%                 colorDict.darkgreen; colorDict.blue; colorDict.cyan]; 
hold on
errorbar(APaxis, average_fittedRate_r1(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(6,:)) % green
errorbar(APaxis, average_fittedRate_r1_mid(:,NC), SEM_fittedRate_r1_mid(:,NC),'Color',ColorChoice(9,:)) % blue
errorbar(APaxis, average_fittedRate_r1_close(:,NC), SEM_fittedRate_r1_close(:,NC),'Color',ColorChoice(10,:)) % dark green

% errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'--','Color',ColorChoice(1,:)) % magenta
% errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'--','Color',ColorChoice(4,:)) % red

ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis - 1 Runt binding site')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('far','middle','close')
%legend('far','middle','close','r0','r3')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_averaged_Different_Position' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_averaged_Different_Position' , '.pdf']);

%% Comparing the 2 binding site constructs
% ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
%                 colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
%                 colorDict.darkgreen; colorDict.blue; colorDict.cyan]; 

hold on
errorbar(APaxis, average_fittedRate_r2(:,NC), SEM_fittedRate_r2(:,NC),'Color',ColorChoice(3,:))
errorbar(APaxis, average_fittedRate_r2_far(:,NC), SEM_fittedRate_r2_far(:,NC),'Color',ColorChoice(7,:))
errorbar(APaxis, average_fittedRate_r2_close(:,NC), SEM_fittedRate_r2_close(:,NC),'Color',ColorChoice(5,:))

% Guide Line
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'--','Color',ColorChoice(1,:)) % magenta
errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'--','Color',ColorChoice(4,:)) % red


ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis - 2 Runt binding sites')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
%legend('2+3','1+3','1+2')
legend('2+3','1+3','1+2','r0','r3')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_r2_averaged_Different_Position_guideline' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_r2_averaged_Different_Position_guideline' , '.pdf']);

%% Plotting the averaged initial RNAP loading rates with SEM
% for all different constructs (mixed sex

hold on
% 0 (No) binding site
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'Color',ColorChoice(1,:)) % magenta
% 3 binding sites
errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'Color',ColorChoice(4,:)) % red
% 1 bining site
errorbar(APaxis, average_fittedRate_r1(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(6,:)) % green
errorbar(APaxis, average_fittedRate_r1_mid(:,NC), SEM_fittedRate_r1_mid(:,NC),'Color',ColorChoice(9,:)) % blue
errorbar(APaxis, average_fittedRate_r1_close(:,NC), SEM_fittedRate_r1_close(:,NC),'Color',ColorChoice(10,:)) % dark green
% % 2 binding sites
errorbar(APaxis, average_fittedRate_r2(:,NC), SEM_fittedRate_r2(:,NC),'Color',ColorChoice(3,:))
errorbar(APaxis, average_fittedRate_r2_far(:,NC), SEM_fittedRate_r2_far(:,NC),'Color',ColorChoice(7,:))
errorbar(APaxis, average_fittedRate_r2_close(:,NC), SEM_fittedRate_r2_close(:,NC),'Color',ColorChoice(5,:))

ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
legend('0','3','1(far)','1(mid)','1(close)',...
        '2(2+3)','2(1+3)','2(1+2)')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_averaged_AllConstruct_r0r3_r1r2variants' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_averaged_AllConstruct_r0r3_r1r2variants' , '.pdf']);

%% Comparing the "1+1 = 2" binding site constructs
% ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
%                 colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
%                 colorDict.darkgreen; colorDict.blue; colorDict.cyan]; 

hold on
errorbar(APaxis, average_fittedRate_r1(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(6,:))
errorbar(APaxis, average_fittedRate_r1_close(:,NC), SEM_fittedRate_r1_close(:,NC),'Color',ColorChoice(10,:))
errorbar(APaxis, average_fittedRate_r2_far(:,NC), SEM_fittedRate_r2_far(:,NC),'Color',ColorChoice(7,:))

% Guide Line
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'--','Color',ColorChoice(1,:)) % magenta
errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'--','Color',ColorChoice(4,:)) % red


ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis - 1 + 1 = 2 (?)')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
%legend('2+3','1+3','1+2')
legend('r1-far','r1-close','1+3','r0','r3')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_close_far_r2_far_averaged_Different_Position_guideline' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_close_far_r2_far_averaged_Different_Position_guideline' , '.pdf']);


%% Comparing the "1+1 = 2" binding site constructs
% ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
%                 colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
%                 colorDict.darkgreen; colorDict.blue; colorDict.cyan]; 

hold on
errorbar(APaxis, average_fittedRate_r1_mid(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(6,:))
errorbar(APaxis, average_fittedRate_r1_close(:,NC), SEM_fittedRate_r1_close(:,NC),'Color',ColorChoice(10,:))
errorbar(APaxis, average_fittedRate_r2(:,NC), SEM_fittedRate_r2(:,NC),'Color',ColorChoice(3,:))
% Guide Line
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'--','Color',ColorChoice(1,:)) % magenta
errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'--','Color',ColorChoice(4,:)) % red


ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis - 1 + 1 = 2 (?)')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
%legend('2+3','1+3','1+2')
legend('r1-mid','r1-close','2+3','r0','r3')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_mid_close_r2_2_+_3_averaged_Different_Position_guideline' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_mid_close_r2_2_+_3_averaged_Different_Position_guideline' , '.pdf']);

%% Comparing the "1+1 = 2" binding site constructs
% ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
%                 colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
%                 colorDict.darkgreen; colorDict.blue; colorDict.cyan]; 

hold on
errorbar(APaxis, average_fittedRate_r1_mid(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(9,:))
errorbar(APaxis, average_fittedRate_r1(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(6,:))
errorbar(APaxis, average_fittedRate_r2_close(:,NC), SEM_fittedRate_r2_close(:,NC),'Color',ColorChoice(10,:))
% Guide Line
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'--','Color',ColorChoice(1,:)) % magenta
errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'--','Color',ColorChoice(4,:)) % red


ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis - 1 + 1 = 2 (?)')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
%legend('2+3','1+3','1+2')
legend('r1-mid','r1-far','1+2','r0','r3')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_mid_far_r2_1_+_2_averaged_Different_Position_guideline' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_r1_mid_far_r2_1_+_2_averaged_Different_Position_guideline' , '.pdf']);

%% Comparing the r3' and r3, also r0, r1
% ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow;...
%                 colorDict.red; colorDict.brown; colorDict.green; colorDict.purple ;...
%                 colorDict.darkgreen; colorDict.blue; colorDict.cyan]; 

hold on
errorbar(APaxis, average_fittedRate_r0(:,NC), SEM_fittedRate_r0(:,NC),'Color',ColorChoice(1,:)) % magenta
errorbar(APaxis, average_fittedRate_r1(:,NC), SEM_fittedRate_r1(:,NC),'Color',ColorChoice(2,:))

errorbar(APaxis, average_fittedRate_r3(:,NC), SEM_fittedRate_r3(:,NC),'Color',ColorChoice(4,:)) % red
errorbar(APaxis, average_fittedRate_r3prime(:,NC), SEM_fittedRate_r3prime(:,NC),'Color',ColorChoice(11,:))


ylim([0 300])
xlim([0.15 0.6])
title('Initial RNAP loading rate over AP axis')
xlabel('AP axis (EL)')
ylabel('Initial RNAP loading rate (AU)')
%legend('2+3','1+3','1+2')
legend('r0','r1','r3','r3 (mutated)')

StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'InitialSlope_r0r1r3_mutated' , '.tif']);
saveas(gcf,[FigPath,filesep,'InitialSlope_r0r1r3_mutated' , '.pdf']);