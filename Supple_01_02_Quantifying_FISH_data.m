function Supple_01_02_Quantifying_FISH_data
%% DESCRIPTION
% We utilize the FISH data from Jeehae Park, Depace lab.
% What we're wondering is whether the normalization that they did was fair.

% Check whether the hkb co-staining + normalization is done or not.
%% Define the directory
%gives the location of these 5 quantities.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
%% Load the P2P-lacZ FISH datasets
FISHPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\hbP2P-LacZ-FISH_LineTraces_Jeehae\wtHbP2-LacZ in situ traces for Yangjoon Kim';
D = dir(FISHPath);
for i=3:length(D)
    FISH_data(i-2,:) = load(['E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\hbP2P-LacZ-FISH_LineTraces_Jeehae\wtHbP2-LacZ in situ traces for Yangjoon Kim\',D(i).name]);
end

%% Plot the data over AP
hold on
for i=1:length(FISH_data(:,1))
    plot(0.005:0.01:0.995, FISH_data(i,:))
    pause
end

title ('FISH pattern over AP')
xlabel('AP axis (EL)')
ylabel('Normalized intensity (AU)')
StandardFigure(gcf,gca)

% Save the figure
%saveas(gcf,'E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\hbP2P-LacZ-FISH_LineTraces_Jeehae\FISH_MS2_Compare_plots\FISH_mRNA_overAP_Jeehae.pdf')

%% Re-normalizing with hkb signal
% I think the FISH intensity profile that I got from Jeehae is normalized
% with maximum intensity, not with hkb signal.
% I'll try to re-normalize the FISH data based on hkb signal (assuming that
% the hkb signal is consistent across embryos, as in Wunderlich and Depace, et al,
% 2014 paper.

% First, for each embryo, get the hkb signal from the posterior 10% of the
% embryo AP axis. (either get mean or maximum)

for embryo=1:length(FISH_data(:,1))
    %hkb signal (mean)
    hkb_mean(embryo) = nanmean(FISH_data(embryo,91:100));
    hkb_max(embryo) = max(FISH_data(embryo,91:100));
    
    FISH_data_renormalized_hkb_Mean(embryo,:) = FISH_data(embryo,:)./hkb_mean(embryo);
    FISH_data_renormalized_hkb_Max(embryo,:) = FISH_data(embryo,:)./hkb_max(embryo);
    
end

%% Plot for checking
hold on
for i=1:length(FISH_data(:,1))
    plot(0.005:0.01:0.995, FISH_data_renormalized_hkb_Mean(i,:))
    pause
end

title ('FISH pattern over AP')
xlabel('AP axis (EL)')
ylabel('Normalized intensity (AU)')
StandardFigure(gcf,gca)

% Save plot
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\FISH_MS2_Compare_plots'
saveas(gcf, [FigPath,filesep,'FISH_P2P_renormalized_hkb_mean.tif'])
saveas(gcf, [FigPath,filesep,'FISH_P2P_renormalized_hkb_mean.pdf'])

%% Averaging the renormalized intensity profile
nEmbryos = length(FISH_data(:,1));
Averaged_FISH_renorm_hkb_Mean = nanmean(FISH_data_renormalized_hkb_Mean);
SEM_FISH_renorm_hkb_Mean = nanstd(FISH_data_renormalized_hkb_Mean,[],1)./sqrt(nEmbryos);

% errorbar(0.005:0.01:0.995, Averaged_FISH_renorm_hkb_Mean, SEM_FISH_renorm_hkb_Mean)

scale_FISH = max(Averaged_FISH_renorm_hkb_Mean);

errorbar(0.005:0.01:0.995,...
         Averaged_FISH_renorm_hkb_Mean./scale_FISH,...
         SEM_FISH_renorm_hkb_Mean./scale_FISH)

title('Renormalized FISH intensity over AP')
xlabel('AP axis (EL)')
ylabel('Renormalized FISH intensity')

StandardFigure(gcf,gca)
    
saveas(gcf, [FigPath,filesep,'FISH_P2P_renormalized_hkb_mean_normalized.tif'])
saveas(gcf, [FigPath,filesep,'FISH_P2P_renormalized_hkb_mean_normalized.pdf'])

%% Compare the renormalized profile vs Jeehae's profile (both averaged over embryos)

hold on
errorbar(0.005:0.01:0.995, nanmean(FISH_data,1), nanstd(FISH_data,[],1))

errorbar(0.005:0.01:0.995,...
         Averaged_FISH_renorm_hkb_Mean./scale_FISH,...
         SEM_FISH_renorm_hkb_Mean./scale_FISH)

ylim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1])
title('Comparison of FISH intensity (different normalization)')
xlabel('AP axis (EL)')
ylabel('Renormalized FISH intensity')
legend('Jeehae','hkb renorm')

StandardFigure(gcf,gca)
    
saveas(gcf, [FigPath,filesep,'FISH_P2P_comparison_DifferentNormalization.tif'])
saveas(gcf, [FigPath,filesep,'FISH_P2P_comparison_DifferentNormalization.pdf'])
%% Save the useful fields
%% Now, let's compare with MS2-MCP data
%% Load the P2P-MS2.V5-lacZ datasets
P2P_MS2V5_data = LoadMS2Sets('P2P-MS2v5-lacZ-36uW')
%P2P_MS2V1_data = LoadMS2Sets('P2P-MS2-lacZ-36uW')
%% Calculate the total mRNA (accumulated mRNA) for P2P-MS2V5 datasets using IntegratemRNA.m
% [TotalProd,TotalProdError,TotalProdN,...
%     MeanTotalProd,SDTotalProd,SETotalProd]=IntegratemRNA(Data,MinParticles,MinEmbryos,varargin)
[TotalProd,TotalProdError,TotalProdN,...
    MeanTotalProd,SDTotalProd,SETotalProd]=IntegratemRNA(P2P_MS2V5_data,3,2);

%% Plot the Mean total Prod (this is total accumulated mRNA over specific NC, divided by the total number of nuclei (ON + Off)
% Since these MeanTotalProd is divided with total # of nuclei, I will
% multiply with the scaling factor to account for the different numbers of
% nuclei in different cycles.

% Scale factor : Assuming that all nuclei undergo doubling at each cycle.
%scale_NC11 = 0.125;
scale_NC12 = 0.25;
scale_NC13 = 0.5;
scale_NC14 = 1;

hold on
%errorbar(0:0.025:1, MeanTotalProd(:,11)*scale_NC11, SETotalProd(:,11)*scale_NC11)
errorbar(0:0.025:1, MeanTotalProd(:,12)*scale_NC12, SETotalProd(:,12)*scale_NC12)
errorbar(0:0.025:1, MeanTotalProd(:,13)*scale_NC13, SETotalProd(:,13)*scale_NC13)
errorbar(0:0.025:1, MeanTotalProd(:,14)*scale_NC14, SETotalProd(:,14)*scale_NC14 )

title('Accumulated mRNA over AP')
xlabel('AP (EL)')
ylabel('Accumulated mRNA (AU)')
legend('NC12','NC13','NC14')

StandardFigure(gcf,gca)

% Save the plot
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\FISH_MS2_Compare_plots'
saveas(gcf, [FigPath,filesep,'AccumulatedmRNA_P2P_MS2MCP_Different_NCs.tif'])
saveas(gcf, [FigPath,filesep,'AccumulatedmRNA_P2P_MS2MCP_Different_NCs.pdf'])
%% Try to calculate the total mRNA over NC12 and NC13

% Scale factor : Assuming that all nuclei undergo doubling at each cycle.
%scale_NC11 = 0.125;
scale_NC12 = 0.5;
scale_NC13 = 1;
scale_NC14 = 0;

MeanTotalProd_NC12NC13 = MeanTotalProd(:,12)*scale_NC12 +...
                            MeanTotalProd(:,13)*scale_NC13; %+ MeanTotalProd(:,14) * scale_NC14; %MeanTotalProd(:,11)*scale_NC11 + 
SETotalProd_NC12NC13 = sqrt (SETotalProd(:,12).^2*scale_NC12 + ...
                                SETotalProd(:,13).^2*scale_NC13); % + SETotalProd(:,14)* scale_NC14);%SETotalProd(:,11).^2*scale_NC11 + 

% Subtract the background (I should do this better)
MeanTotalProd_NC12NC13 = MeanTotalProd_NC12NC13 - min(MeanTotalProd_NC12NC13);
                            
MeanTotalProd_NC12NC13_Normalized = MeanTotalProd_NC12NC13 ./ max(MeanTotalProd_NC12NC13);
SETotalProd_NC12NC13_Normalized = SETotalProd_NC12NC13./max(MeanTotalProd_NC12NC13);

%% plot TotalProd from individual embryos, as well as the mean (for nc12 + nc13, for now)
APaxis = 0:0.025:1;

% Scale factor : Assuming that all nuclei undergo doubling at each cycle.
%scale_NC11 = 0.125;
scale_NC12 = 0.5;
scale_NC13 = 1;
scale_NC14 = 0;

TotalmRNA_individual = nan(length(TotalProd(:,1,1)), 41);
TotalmRNA_individual_Error = nan(length(TotalProd(:,1,1)), 41);

% Calculate the TotalmRNA and error for nc12 and nc13
for i=1:length(TotalProd(:,1,1))
    TotalmRNA_individual(i,:) = TotalProd(i,:,12)*scale_NC12 + TotalProd(i,:,13)*scale_NC13;
    TotalmRNA_individual_Error(i,:) = sqrt(TotalProdError(i,:,12).^2*scale_NC12 + ...
                           TotalProdError(i,:,13).^2*scale_NC13);
end

hold on
for i=1:length(TotalProd(:,1,1))
    h(i) = errorbar(APaxis, TotalmRNA_individual(i,:),  TotalmRNA_individual_Error(i,:))
    pause
end

h(length(TotalProd(:,1,1))+1) = errorbar(APaxis, MeanTotalProd_NC12NC13, SETotalProd_NC12NC13)

title('Accumulated mRNA (nc12+nc13) : P2P-MS2.V5')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend(h(length(TotalProd(:,1,1))+1),'Mean')
StandardFigure(gcf,gca)

% Save the figure
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\FISH_MS2_Compare_plots'
% saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_MS2V5_Individual_Mean.tif'])
% saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_MS2V5_Individual_Mean.pdf'])
%% Normalize the Mean Total Prod for NC13 and NC14
MeanTotalProd_NC13_Normalized = MeanTotalProd(:,13) ./ max(MeanTotalProd(:,13));
SETotalProd_NC13_Normalized = SETotalProd(:,13)./ max(MeanTotalProd(:,13));

MeanTotalProd_NC14_Normalized = MeanTotalProd(:,14) ./ max(MeanTotalProd(:,14));
SETotalProd_NC14_Normalized = SETotalProd(:,14) ./ max(MeanTotalProd(:,14));


%% Compare the MS2-MCP and FISH

%% Load the in situ data
load ('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon\ProcessedProfile_averaged.mat')
% This is processed by Process_InSitu_data.m and Supple_02.m 
% Averaged_Intensity_rN, SEM_Intensity_rN are the fields.

%% Plot
% FISH (re-normalized with hkb signal)
hold on
errorbar(0.005:0.01:0.995,...
         Averaged_FISH_renorm_hkb_Mean./scale_FISH,...
         SEM_FISH_renorm_hkb_Mean./scale_FISH)
     
% In Situ
% I'll just use the hb P2-evePr (r0) for now, since we don't have any in
% situ data from the P2P. Maybe we can grab the smFISH data of P2P?
errorbar(0.005:0.01:0.995,...
         Averaged_Intensity_r0./max(Averaged_Intensity_r0),...
         SEM_Intensity_r0./max(Averaged_Intensity_r0))

ylim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1])
title('Comparison of FISH intensity (different normalization)')
xlabel('AP axis (EL)')
ylabel('Renormalized FISH intensity')
legend('FISH-hkb renorm (P2P)','in situ(P2-evePr)')

StandardFigure(gcf,gca)
% Save plots

%% Compare in situ, FISH, and MS2
% FISH (re-normalized with hkb signal)
hold on
errorbar(0.005:0.01:0.995,...
         Averaged_FISH_renorm_hkb_Mean./scale_FISH,...
         SEM_FISH_renorm_hkb_Mean./scale_FISH)
     
% In Situ
% I'll just use the hb P2-evePr (r0) for now, since we don't have any in
% situ data from the P2P. Maybe we can grab the smFISH data of P2P?
errorbar(0.005:0.01:0.995,...
         Averaged_Intensity_r0./max(Averaged_Intensity_r0),...
         SEM_Intensity_r0./max(Averaged_Intensity_r0))

% MS2-MCP
APaxis = 0.0125:0.025:1.0125 % For a better registration. The last AP bin is practically nothing.
errorbar(APaxis,...
        MeanTotalProd_NC12NC13_Normalized*1.4,...
        SETotalProd_NC12NC14_Normalized)

ylim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1])
title('Accumulated mRNA')
xlabel('AP axis (EL)')
ylabel('Normalized intensity')
legend('FISH (P2P)','in situ(P2-evePr)', 'MS2-MCP(P2P)')

StandardFigure(gcf,gca)


% Save the plot
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\FISH_MS2_Compare_plots'
saveas(gcf, [FigPath,filesep,'AccumulatedmRNA_Comparison of All Methods.tif'])
saveas(gcf, [FigPath,filesep,'AccumulatedmRNA_Comparison of All Methods.pdf'])

%% Save the useful fields
%% FISH (after re-normalization)
% Save FISH data
savedVariables_FISH = {};
savedVariables_FISH = [savedVariables_FISH,...
                                    'Averaged_FISH_renorm_hkb_Mean',...
                                    'SEM_FISH_renorm_hkb_Mean',...
                                    'FISH_data',...
                                    'hkb_mean'];
                                
% Save the variables into .mat file.
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
save([FilePath,filesep,'FISH_P2P_Jeehae.mat'],...
    savedVariables_FISH{:},'-v7.3');

%% MS2-MCP (P2P)
% Save the MS2-MCP (P2P-MS2V5) data
savedVariables_MS2 = {};
savedVariables_MS2 = [savedVariables_MS2,...
                                    'MeanTotalProd_NC12NC13_Normalized',...
                                    'SETotalProd_NC12NC13_Normalized',...
                                    'MeanTotalProd'];
                                
% Save the variables into .mat file.
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
save([FilePath,filesep,'P2P-MS2V5-Liz.mat'],...
    savedVariables_MS2{:},'-v7.3');

%% MS2-MCP : for r0,1,2,3 (late NC14)
% Save the MS2-MCP (P2P-MS2V5) data
savedVariables_MS2MCP = {};
savedVariables_MS2MCP = [savedVariables_MS2MCP,...
                            'Averaged_integratedmRNA_r0', 'SEM_integratedmRNA_r0'...
                            'Averaged_integratedmRNA_r1', 'SEM_integratedmRNA_r1'...
                            'Averaged_integratedmRNA_r2', 'SEM_integratedmRNA_r2'...
                            'Averaged_integratedmRNA_r3', 'SEM_integratedmRNA_r3'...
                            'Averaged_integratedmRNA_r3prime', 'SEM_integratedmRNA_r3prime'];

% Save the variables into .mat file.
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';
save([FilePath,filesep,'r0123_MS2V5-MCP_lateNC14.mat'],...
    savedVariables_MS2MCP{:},'-v7.3');

end