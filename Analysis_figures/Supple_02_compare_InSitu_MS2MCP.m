function Supple_02_compare_InSitu_MS2MCP
% DESCRIPTION
% This script is for comparison of In Situ hybridization and MS2-MCP 
% in terms of accumulated cytoplasmic mRNA pattern. 
% As the first step, we will compare Steve Small lab's data for P2-r0,1,2,3, with our
% P2(r0,1,2,3)-MS2-MCP data.

%% Load the In Situ datasets
InSituPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon';
r0Path = [InSituPath,filesep,'hbP2\NC14'];
r1Path = [InSituPath,filesep,'hbP2+1Run\NC14'];
r2Path = [InSituPath,filesep,'hbP2+2Run\NC14'];
r3Path = [InSituPath,filesep,'hbP2+3Run\NC14'];

% Optional
r0Path_earlyNC14 = [InSituPath,filesep,'hbP2\earlyNC14'];
%% Load specific embryos
% Use dir function to read out all names in the directory.
D = dir(r0Path_earlyNC14);
DataPath = r0Path_earlyNC14;

%% Loop through all images (embryos)
%k=1; % Count the number of imagesath;

AP = 0.005:0.01:0.995;
    
for i=3:length(D)
    % import the raw image
    EmbryoName = D(i).name;
    InSituEmbryo_raw = imread([DataPath,filesep, EmbryoName]);
    % Get the raw, binned, smoothened, and backgroud-subtracted intensity
    % profile over the defined AP axis.
    [Intensity_raw,...
            nBins, Intensity_mean,  Intensity_std, ...
            window_smooth, Intensity_mean_smooth,...
            background_intensity, Intensity_mean_smooth_BGsubtracted, X, Y] =...
                        Process_InSitu_data(InSituEmbryo_raw, EmbryoName, DataPath);
        
    % Save all the processed data into one structure, this is similar to making CompiledParticles.mat    
    ProcessedInSituData(i-2).EmbryoName = EmbryoName;
    
    % Raw intensity profile over AP axis line.        
    ProcessedInSituData(i-2).Intensity_raw = Intensity_raw;
    
    % Binning
    ProcessedInSituData(i-2).nBins = nBins;
    ProcessedInSituData(i-2).Intensity_mean = Intensity_mean;
    ProcessedInSituData(i-2).Intensity_std = Intensity_std;
    
    % Smoothing
    ProcessedInSituData(i-2).window_smooth = window_smooth;
    ProcessedInSituData(i-2).Intensity_mean_smooth = Intensity_mean_smooth;
    
    % Background-subtraction
    ProcessedInSituData(i-2).Background = background_intensity;
    ProcessedInSituData(i-2).Intensity_mean_smooth_BGsubtracted = Intensity_mean_smooth_BGsubtracted;
    
    % Head and Tail points
    ProcessedInSituData(i-2).X = X;
    ProcessedInSituData(i-2).Y = Y;
    %k=k+1;
end

%% Save the fields
% I want to save the processed data
% 1) the raw embryo image with AP axis defined (DONE)
% 2) Intensity (raw, binned, smoothened)

save([DataPath,filesep,'ProcessedInSituData',DataPath(end-8:end-5),'.mat'],'ProcessedInSituData')

%% Get the boundary features from processed data

% Load the AP-defined datasets

AP = 0.005:0.01:0.995;

for i=1:length(ProcessedInSituData)
    clf
    % Define the fields from the ProcessedInSituData
    EmbryoName = ProcessedInSituData(i).EmbryoName
    
    Intensity_raw = ProcessedInSituData(i).Intensity_raw;
    Intensity_mean = ProcessedInSituData(i).Intensity_mean;
    Intensity_std = ProcessedInSituData(i).Intensity_std;
    Intensity_mean_smooth = ProcessedInSituData(i).Intensity_mean_smooth;
    Background = ProcessedInSituData(i).Background;
    Intensity_mean_smooth_BGsubtracted = ProcessedInSituData(i).Intensity_mean_smooth_BGsubtracted;
    
    % 1) Get the boundary features, as our way, half-maximum
    % First, let's use our standard, which is getting the half-maximum, and
    % tangential line at that point.
    % [Max, Min, DynamicRange, BoundaryPosition, Slope] =...
    %                         getGradientFeatures (AP,Intensity_mean_smooth);                    
    [Max, Min, DynamicRange, BoundaryPosition, Slope] =...
                            getGradientFeatures (AP,Intensity_mean_smooth_BGsubtracted);
                        
    [BoundaryPosition_inflection, Slope_inflection] = ...
                        Get_Boundary_inflection_points(AP,Intensity_mean_smooth_BGsubtracted)
    
    % Boundary features
    ProcessedInSituData(i).BoundaryPosition = BoundaryPosition;
    ProcessedInSituData(i).Slope = Slope;
    
    ProcessedInSituData(i).BoundaryPosition_inflection = BoundaryPosition_inflection;
    ProcessedInSituData(i).Slope_inflection = Slope_inflection;
    



                    
    % Plot the boundary features on top of the profile (half-maximum)
    hold on
    % Raw intensity
    plot((1:length(Intensity_raw))/length(Intensity_raw) ,Intensity_raw)
    % Averaged pixel intensity, std (over AP)- binned
    errorbar(0.005:0.01:0.995,Intensity_mean,Intensity_std)
    % Smoothened over window_smooth
    plot(0.005:0.01:0.995,Intensity_mean_smooth)
    
    plot(BoundaryPosition, Intensity_mean_smooth(AP == BoundaryPosition),'o','MarkerSize',10)
    plot(BoundaryPosition_inflection,Intensity_mean_smooth(AP == BoundaryPosition_inflection),'o','MarkerSize',10)
    
    % Bar line for the Background
    plot(AP, ones(size(AP))*Background,'k')
    
    % tangential line
    plot(AP,Slope*(AP-BoundaryPosition) + Intensity_mean_smooth(AP == BoundaryPosition))
    

    %h(7) = plot(AP,Slope*(AP-BoundaryPosition_inflection) + Intensity_mean_smooth(AP == BoundaryPosition_inflection))
    
    ylim([0 max(Intensity_mean_smooth)+50])
    title('Intensity from In Situ embryo')
    xlabel('AP (EL)')
    ylabel('Intensity (AU)')
    legend('Raw','Binned','Smoothened','half-maximum','Inflection')
    StandardFigure(gcf,gca)

    % Save the result figure
    saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'.tif'])
    saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'.pdf'])
end

%% Compare the In Situ intensity (background subtracted) of multiple embryos

AP = 0.005:0.01:0.995;

hold on
for i=1:length(ProcessedInSituData)
    Intensity_mean_smooth_BGsubtracted(i,:) = ProcessedInSituData(i).Intensity_mean_smooth_BGsubtracted;
    plot(AP, ProcessedInSituData(i).Intensity_mean_smooth_BGsubtracted)
    pause
end

title('Individual In Situ intensity profile over AP')
xlabel('AP axis (EL)')
ylabel('In Situ intensity (AU)')
StandardFigure(gcf,gca)

% Save the result figure
saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,DataPath(end-8:end-5),'AllEmbryos_BGsubtracted.tif'])
saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,DataPath(end-8:end-5),'AllEmbryos_BGsubtracted.pdf'])

%save([DataPath,filesep,'ProcessedInSituData',DataPath(end-8:end-5),'.mat'],'ProcessedInSituData')
%% Average the BG subtracted Intensity profile (of all embryos)
Averaged_Intensity_BGsubtracted = nanmean(Intensity_mean_smooth_BGsubtracted);
SEM_Intensity_BGsubtracted = nanstd(Intensity_mean_smooth_BGsubtracted)./sqrt(length(ProcessedInSituData));
errorbar(AP, Averaged_Intensity_BGsubtracted, SEM_Intensity_BGsubtracted)

title('Averaged In Situ intensity profile over AP')
xlabel('AP axis (EL)')
ylabel('In Situ intensity (AU)')
StandardFigure(gcf,gca)

% Save the result figure
saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,DataPath(end-8:end-5),'_averaged.tif'])
saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,DataPath(end-8:end-5),'_averaged.pdf'])

%% Save the averaged fields
Averaged_insitu_profile.Averaged_Intensity_BGsubtracted = Averaged_Intensity_BGsubtracted;
Averaged_insitu_profile.SEM_Intensity_BGsubtracted = SEM_Intensity_BGsubtracted;

save([DataPath,filesep,'Averaged_InSitu_profile','_earlyNC14','.mat'],'Averaged_insitu_profile')

%% Save the values for further comparison with other constructs
Averaged_Intensity_r3 = Averaged_Intensity_BGsubtracted;
SEM_Intensity_r3 = SEM_Intensity_BGsubtracted;

load ('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon\ProcessedProfile_averaged.mat')
save('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon\ProcessedProfile_averaged.mat',...
            'Averaged_Intensity_r0','SEM_Intensity_r0',...
            'Averaged_Intensity_r1','SEM_Intensity_r1',...
            'Averaged_Intensity_r2','SEM_Intensity_r2',...
            'Averaged_Intensity_r3','SEM_Intensity_r3')

%% Part2. Start from here in case I already processed the raw data
%% Load the datasets
load ('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon\ProcessedProfile_averaged.mat')
%% Plot the averaged In Situ intensity profile for r0,1,2,3

AP = 0:0.01:0.99;
AP = 0.005:0.01:0.995;

hold on
errorbar(AP, Averaged_Intensity_r0, SEM_Intensity_r0,'color',ColorChoice(1,:))
errorbar(AP, Averaged_Intensity_r1, SEM_Intensity_r1,'color',ColorChoice(2,:))
errorbar(AP, Averaged_Intensity_r2, SEM_Intensity_r2,'color',ColorChoice(3,:))
errorbar(AP, Averaged_Intensity_r3, SEM_Intensity_r3,'color',ColorChoice(4,:))

ylim([0 160])
xticks([0 0.2 0.4 0.6 0.8 1])

title('{\it in situ} intensity')
xlabel('AP axis (EL)')
ylabel('{\it in situ} intensity (AU)')
legend('r0','r1','r2','r3')
StandardFigure(gcf,gca)

% % Save the result figure
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\OpposingGradientsFigures\Data\InSitu_MS2_compare_plots\ProcessedResults_InSitu'
saveas(gcf,[figPath, filesep, 'r0123_averaged_intensity_InSitu.tif'])
saveas(gcf,[figPath, filesep, 'r0123_averaged_intensity_InSitu.pdf'])

%% Plot the averaged In Situ intensity profile for r0,1,2,3 (normalized)

hold on
errorbar(AP, Averaged_Intensity_r0./max(Averaged_Intensity_r0), SEM_Intensity_r0./max(Averaged_Intensity_r0))
errorbar(AP, Averaged_Intensity_r1./max(Averaged_Intensity_r1), SEM_Intensity_r1./max(Averaged_Intensity_r1))
errorbar(AP, Averaged_Intensity_r2./max(Averaged_Intensity_r2), SEM_Intensity_r2./max(Averaged_Intensity_r2))
errorbar(AP, Averaged_Intensity_r3./max(Averaged_Intensity_r3), SEM_Intensity_r3./max(Averaged_Intensity_r3))

title('In Situ intensity (averaged) over AP')
xlabel('AP axis (EL)')
ylabel('In Situ intensity (AU)')
legend('r0','r1','r2','r3')
StandardFigure(gcf,gca)

% Save the result figure
saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon\r0123_averaged_intensity_InSitu_Normalized.tif'])
saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon\r0123_averaged_intensity_InSitu_Normalized.pdf'])

%% %%%%%%%%%%%%%%%%%% Let's compare this with MS2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AccumulatedmRNA_mixed_sex';
% % First, use Average_TotalmRNAProd function (custom-written). 
% % Put some more description later on how it works.
% % This 'Average_TotalmRNAProd' uses the mat file produced by
% % 'AverageDatasets.m'. 
% [TotalmRNA_averaged_r0, TotalmRNA_SEM_r0] = Average_TotalmRNAProd('r0-', DataPath);
% [TotalmRNA_averaged_r1, TotalmRNA_SEM_r1] = Average_TotalmRNAProd('r1-', DataPath);
% [TotalmRNA_averaged_r2, TotalmRNA_SEM_r2] = Average_TotalmRNAProd('r2-', DataPath);
% [TotalmRNA_averaged_r3, TotalmRNA_SEM_r3] = Average_TotalmRNAProd('r3-', DataPath);
% 
% % save the fields
% %% Plot AccumulatedmRNA (MS2) altogether
% APaxis = 0:0.025:1;
% hold on
% errorbar(APaxis, TotalmRNA_averaged_r0, TotalmRNA_SEM_r0)
% errorbar(APaxis, TotalmRNA_averaged_r1, TotalmRNA_SEM_r1)
% errorbar(APaxis, TotalmRNA_averaged_r2, TotalmRNA_SEM_r2)
% errorbar(APaxis, TotalmRNA_averaged_r3, TotalmRNA_SEM_r3)
% 
% title('Accumulated mRNA (MS2, averaged) over AP')
% xlabel('AP axis (EL)')
% ylabel('Accumulated mRNA (AU)')
% legend('r0','r1','r2','r3')
% StandardFigure(gcf,gca)

% Save the result figure
% saveas(gcf,[DataPath,filesep,'AccumulatedmRNA_r0123_MS2_averaged','.tif'])
% saveas(gcf,[DataPath,filesep,'AccumulatedmRNA_r0123_MS2_averaged','.pdf'])

%% Optional (Calculate the totalmRNA using IntegratemRNA.m script)
% This has benefit that we can control how many embryos, particles we can
% use as a threshold to cut out.
r0Data = LoadMS2Sets('r0');
r1Data = LoadMS2Sets('r1');
r2Data = LoadMS2Sets('r2');
r3Data = LoadMS2Sets('r3');
r3primeData = LoadMS2Sets('r3prime');

%% Calculate the TotalmRNA using IntegratemRNA script
[TotalProd_r0,TotalProdError_r0,TotalProdN_r0,...
    MeanTotalProd_r0,SDTotalProd_r0,SETotalProd_r0]=IntegratemRNA(r0Data,2,2)

[TotalProd_r1,TotalProdError_r1,TotalProdN_r1,...
    MeanTotalProd_r1,SDTotalProd_r1,SETotalProd_r1]=IntegratemRNA(r1Data,2,2)

[TotalProd_r2,TotalProdError_r2,TotalProdN_r2,...
    MeanTotalProd_r2,SDTotalProd_r2,SETotalProd_r2]=IntegratemRNA(r2Data,2,2)

[TotalProd_r3,TotalProdError_r3,TotalProdN_r3,...
    MeanTotalProd_r3,SDTotalProd_r3,SETotalProd_r3]=IntegratemRNA(r3Data,2,2)

[TotalProd_r3prime,TotalProdError_r3prime,TotalProdN_r3prime,...
    MeanTotalProd_r3prime,SDTotalProd_r3prime,SETotalProd_r3prime]=IntegratemRNA(r3primeData,2,2)


%% Plot the IntegratedmRNA from all different constructs

% Define the NC
NC=13;

% r0
figure_integratedmRNA_r0 = figure;
hold on
for i=1:length(TotalProd_r0(:,1,1))
    errorbar(0:0.025:1, TotalProd_r0(i,:,NC), TotalProdError_r0(i,:,NC))
end

% r1
figure_integratedmRNA_r1 = figure;
hold on
for i=1:length(TotalProd_r1(:,1,1))
    errorbar(0:0.025:1, TotalProd_r1(i,:,NC), TotalProdError_r1(i,:,NC))
end

% r2
figure_integratedmRNA_r2 = figure;
hold on
for i=1:length(TotalProd_r2(:,1,1))
    errorbar(0:0.025:1, TotalProd_r2(i,:,NC), TotalProdError_r2(i,:,NC))
end

% r3
figure_integratedmRNA_r3 = figure;
hold on
for i=1:length(TotalProd_r3(:,1,1))
    errorbar(0:0.025:1, TotalProd_r3(i,:,NC), TotalProdError_r3(i,:,NC))
end

%% Calculate the amount during NC13 to NC14
Scale_NC13 = 0.5;
Scale_NC14 = 1;

integratedmRNA_r0 = TotalProd_r0(:,:,13)*Scale_NC13 + ...
                    TotalProd_r0(:,:,14)*Scale_NC14;
                
integratedmRNA_r1 = TotalProd_r1(:,:,13)*Scale_NC13 + ...
                    TotalProd_r1(:,:,14)*Scale_NC14;
                
integratedmRNA_r2 = TotalProd_r2(:,:,13)*Scale_NC13 + ...
                    TotalProd_r2(:,:,14)*Scale_NC14;
                
integratedmRNA_r3 = TotalProd_r3(:,:,13)*Scale_NC13 + ...
                    TotalProd_r3(:,:,14)*Scale_NC14;
                
integratedmRNA_r3prime = TotalProd_r3prime(:,:,13)*Scale_NC13 + ...
                    TotalProd_r3prime(:,:,14)*Scale_NC14;
                
% Error estimation
integratedmRNA_error_r0 = sqrt(TotalProd_r0(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r0(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r1 = sqrt(TotalProd_r1(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r1(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r2 = sqrt(TotalProd_r2(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r2(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r3 = sqrt(TotalProd_r3(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r3(:,:,14).^2*Scale_NC14);
                        
integratedmRNA_error_r3prime = sqrt(TotalProd_r3prime(:,:,13).^2*Scale_NC13 + ...
                            TotalProd_r3prime(:,:,14).^2*Scale_NC14);

%% Plot the IntegratedmRNA from all different constructs (NC13 + NC14)


% r0
figure_integratedmRNA_r0 = figure;
hold on
for i=1:length(integratedmRNA_r0(:,1))
    errorbar(0:0.025:1, integratedmRNA_r0(i,:), integratedmRNA_error_r0(i,:))
end

% r1
figure_integratedmRNA_r1 = figure;
hold on
for i=1:length(integratedmRNA_r1(:,1))
    errorbar(0:0.025:1, integratedmRNA_r1(i,:), integratedmRNA_error_r1(i,:))
end       

% r2
figure_integratedmRNA_r2 = figure;
hold on
for i=1:length(integratedmRNA_r2(:,1))
    errorbar(0:0.025:1, integratedmRNA_r2(i,:), integratedmRNA_error_r2(i,:))
end

% r3
figure_integratedmRNA_r3 = figure;
hold on
for i=1:length(integratedmRNA_r3(:,1))
    errorbar(0:0.025:1, integratedmRNA_r3(i,:), integratedmRNA_error_r3(i,:))
end

% r3 prime
figure_integratedmRNA_r3prime = figure;
hold on
for i=1:length(integratedmRNA_r3prime(:,1))
    errorbar(0:0.025:1, integratedmRNA_r3prime(i,:), integratedmRNA_error_r3prime(i,:))
end

%% Optional 
%% Averaging over multiple embryos
Averaged_integratedmRNA_r0 = nanmean(integratedmRNA_r0);
Averaged_integratedmRNA_r1 = nanmean(integratedmRNA_r1);
Averaged_integratedmRNA_r2 = nanmean(integratedmRNA_r2);
Averaged_integratedmRNA_r3 = nanmean(integratedmRNA_r3);
Averaged_integratedmRNA_r3prime = nanmean(integratedmRNA_r3prime);

SEM_integratedmRNA_r0 = nanstd(integratedmRNA_r0,[],1)./sqrt(length(integratedmRNA_r0(:,1)));
SEM_integratedmRNA_r1 = nanstd(integratedmRNA_r1,[],1)./sqrt(length(integratedmRNA_r1(:,1)));
SEM_integratedmRNA_r2 = nanstd(integratedmRNA_r2,[],1)./sqrt(length(integratedmRNA_r2(:,1)));
SEM_integratedmRNA_r3 = nanstd(integratedmRNA_r3,[],1)./sqrt(length(integratedmRNA_r3(:,1)));
SEM_integratedmRNA_r3prime = nanstd(integratedmRNA_r3prime,[],1)./sqrt(length(integratedmRNA_r3prime(:,1)));

%% Plot to check (averaged accumulated mRNA profile)
APaxis = 0:0.025:1;

hold on
errorbar(APaxis, Averaged_integratedmRNA_r0, SEM_integratedmRNA_r0)
errorbar(APaxis, Averaged_integratedmRNA_r1, SEM_integratedmRNA_r1)
errorbar(APaxis, Averaged_integratedmRNA_r2, SEM_integratedmRNA_r2)
errorbar(APaxis, Averaged_integratedmRNA_r3, SEM_integratedmRNA_r3)
errorbar(APaxis, Averaged_integratedmRNA_r3prime, SEM_integratedmRNA_r3prime)

title('Accumulated mRNA')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('0','1','2','3','3(mutated)')

StandardFigure(gcf,gca)

% Save the plots
% File path
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA'
saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_AllConstructs(SmallLab)-MS2MCP-Averaged.tif'])
saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_AllConstructs(SmallLab)-MS2MCP-Averaged.pdf'])


%% Plot after normalization
scale_r0 = max(Averaged_integratedmRNA_r0);
scale_r1 = max(Averaged_integratedmRNA_r1);
scale_r2 = max(Averaged_integratedmRNA_r2);
scale_r3 = max(Averaged_integratedmRNA_r3);
scale_r3prime = max(Averaged_integratedmRNA_r3prime);

hold on
errorbar(APaxis, Averaged_integratedmRNA_r0./scale_r0, SEM_integratedmRNA_r0./scale_r0)
errorbar(APaxis, Averaged_integratedmRNA_r1./scale_r1, SEM_integratedmRNA_r1./scale_r1)
errorbar(APaxis, Averaged_integratedmRNA_r2./scale_r2, SEM_integratedmRNA_r2./scale_r2)
errorbar(APaxis, Averaged_integratedmRNA_r3./scale_r3, SEM_integratedmRNA_r3./scale_r3)
errorbar(APaxis, Averaged_integratedmRNA_r3prime./scale_r3prime, SEM_integratedmRNA_r3prime./scale_r3prime)

title('Accumulated mRNA')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('0','1','2','3','3(mutated)')

StandardFigure(gcf,gca)

% Save the plots
% File path
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\AccumulatedmRNA'
saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_AllConstructs(SmallLab)-MS2MCP-Normalized.tif'])
saveas(gcf,[FigPath,filesep,'AccumulatedmRNA_AllConstructs(SmallLab)-MS2MCP-Normalized.pdf'])

%% Compare the In Situ vs MS2 profile (Normalized)
% Caveat : I'm using the MS2-integrated mRNA calculated by AverageDatasets.m

AP = 0.005:0.01:0.995;
APaxis = 0:0.025:1;

% r0
figure_compare_r0 = figure(1)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r0./max(Averaged_Intensity_r0), SEM_Intensity_r0./max(Averaged_Intensity_r0))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r0./max(TotalmRNA_averaged_r0), TotalmRNA_SEM_r0./max(TotalmRNA_averaged_r0))

title('Accumulated mRNA over AP - r0')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r1
figure_compare_r1 = figure(2)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r1./max(Averaged_Intensity_r1), SEM_Intensity_r1./max(Averaged_Intensity_r1))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r1./max(TotalmRNA_averaged_r1), TotalmRNA_SEM_r1./max(TotalmRNA_averaged_r1))

title('Accumulated mRNA over AP - r1')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r2
figure_compare_r2 = figure(3)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r2./max(Averaged_Intensity_r2), SEM_Intensity_r2./max(Averaged_Intensity_r2))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r2./max(TotalmRNA_averaged_r2), TotalmRNA_SEM_r2./max(TotalmRNA_averaged_r2))

title('Accumulated mRNA over AP - r2')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r3
figure_compare_r3 = figure(4)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r3./max(Averaged_Intensity_r3), SEM_Intensity_r3./max(Averaged_Intensity_r3))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r3./max(TotalmRNA_averaged_r3), TotalmRNA_SEM_r3./max(TotalmRNA_averaged_r3))

title('Accumulated mRNA over AP - r3')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

%% Save all plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InSitu_MS2_compare_plots\Normalized_profile_comparison';

saveas(figure_compare_r0, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r0.tif'])
saveas(figure_compare_r0, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r0.pdf'])

saveas(figure_compare_r1, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r1.tif'])
saveas(figure_compare_r1, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r1.pdf'])

saveas(figure_compare_r2, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r2.tif'])
saveas(figure_compare_r2, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r2.pdf'])

saveas(figure_compare_r3, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r3.tif'])
saveas(figure_compare_r3, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r3.pdf'])

%% Compare the In Situ vs MS2 profile (Assuming some saturation of In Situ)
%% step1. Deciding the cut-off line
% From the In Situ profile of r0,1,2,3, it seems that r0 has a broad peak
% up to 40% of the embryo, thus I'll make a rough assumption that the 40%
% of the r0 signal is the saturation point.

% This cutoff point can be tuned such that the 
%CutOff = TotalmRNA_averaged_r0(find(APaxis==0.4));
%CutOff = TotalmRNA_averaged_r0(find(APaxis==0.3500));
CutOff = TotalmRNA_averaged_r0(20); % 42.5% of EL

APaxis = 0:0.025:1;
hold on
errorbar(APaxis, TotalmRNA_averaged_r0, TotalmRNA_SEM_r0)
errorbar(APaxis, TotalmRNA_averaged_r1, TotalmRNA_SEM_r1)
errorbar(APaxis, TotalmRNA_averaged_r2, TotalmRNA_SEM_r2)
errorbar(APaxis, TotalmRNA_averaged_r3, TotalmRNA_SEM_r3)

plot(APaxis, ones(size(APaxis))* CutOff)

title('Accumulated mRNA (MS2, averaged) over AP')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (AU)')
legend('r0','r1','r2','r3','Saturation')
StandardFigure(gcf,gca)

% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InSitu_MS2_compare_plots\CutOff_Normalized_profile_comparison';
% saveas(gcf, [FigPath,filesep,'TotalmRNA_MS2_CutOff.tif'])
% saveas(gcf, [FigPath,filesep,'TotalmRNA_MS2_CutOff.pdf'])
%% step2. Cap the amplitude of TotalmRNA using the CutOff value.
% First, define the values to be modified.
TotalmRNA_averaged_r0_CutOff = TotalmRNA_averaged_r0;
TotalmRNA_averaged_r1_CutOff = TotalmRNA_averaged_r1;
TotalmRNA_averaged_r2_CutOff = TotalmRNA_averaged_r2;
TotalmRNA_averaged_r3_CutOff = TotalmRNA_averaged_r3;

TotalmRNA_averaged_r0_CutOff(TotalmRNA_averaged_r0_CutOff >= CutOff) = CutOff;
TotalmRNA_averaged_r1_CutOff(TotalmRNA_averaged_r1_CutOff >= CutOff) = CutOff;
TotalmRNA_averaged_r2_CutOff(TotalmRNA_averaged_r2_CutOff >= CutOff) = CutOff;
TotalmRNA_averaged_r3_CutOff(TotalmRNA_averaged_r3_CutOff >= CutOff) = CutOff;

TotalmRNA_SEM_r0_CutOff = TotalmRNA_SEM_r0;
TotalmRNA_SEM_r1_CutOff = TotalmRNA_SEM_r1;
TotalmRNA_SEM_r2_CutOff = TotalmRNA_SEM_r2;
TotalmRNA_SEM_r3_CutOff = TotalmRNA_SEM_r3;

TotalmRNA_SEM_r0_CutOff(TotalmRNA_averaged_r0_CutOff >= CutOff) = 0;
TotalmRNA_SEM_r1_CutOff(TotalmRNA_averaged_r1_CutOff >= CutOff) = 0;
TotalmRNA_SEM_r2_CutOff(TotalmRNA_averaged_r2_CutOff >= CutOff) = 0;
TotalmRNA_SEM_r3_CutOff(TotalmRNA_averaged_r3_CutOff >= CutOff) = 0;

%% step3. Plotting the re-normalized curves together

AP = 0.005:0.01:0.995;
APaxis = 0:0.025:1;

% r0
figure_compare_r0_CutOff = figure(1)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r0./max(Averaged_Intensity_r0), SEM_Intensity_r0./max(Averaged_Intensity_r0))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r0_CutOff./CutOff, TotalmRNA_SEM_r0_CutOff./CutOff)

title('Accumulated mRNA over AP - r0')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r1
figure_compare_r1_CutOff = figure(2)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r1./max(Averaged_Intensity_r1), SEM_Intensity_r1./max(Averaged_Intensity_r1))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r1_CutOff./CutOff, TotalmRNA_SEM_r1_CutOff./CutOff)

title('Accumulated mRNA over AP - r1')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r2
figure_compare_r2_CutOff = figure(3)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r2./max(Averaged_Intensity_r2), SEM_Intensity_r2./max(Averaged_Intensity_r2))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r2_CutOff./CutOff,TotalmRNA_SEM_r2_CutOff./CutOff)

title('Accumulated mRNA over AP - r2')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r3
figure_compare_r3_CutOff = figure(4)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r3./max(Averaged_Intensity_r3), SEM_Intensity_r3./max(Averaged_Intensity_r3))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r3_CutOff./CutOff, TotalmRNA_SEM_r3_CutOff./CutOff)

title('Accumulated mRNA over AP - r3')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

%% Save plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InSitu_MS2_compare_plots\CutOff_Normalized_profile_comparison';

saveas(figure_compare_r0_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r0.tif'])
saveas(figure_compare_r0_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r0.pdf'])

saveas(figure_compare_r1_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r1.tif'])
saveas(figure_compare_r1_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r1.pdf'])

saveas(figure_compare_r2_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r2.tif'])
saveas(figure_compare_r2_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r2.pdf'])

saveas(figure_compare_r3_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r3.tif'])
saveas(figure_compare_r3_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r3.pdf'])

%% CutOff value : Finding the optimal value for each construct
% Added by YJK on 9/24/2019
%% step1. Find the optimal CutOff value for each construct
% The idea is to move the cutoff bar from 0 to maximum of MS2 to find the
% one that minimizes the squared sum of (in situ - adjusted MS2), basically
% lsqnonlin (least square)
% Caveat : Let's only compare the region of 0.2 to 0.6

% r0
CapRange = 0:100:10000;

CapIntensity_r0 = nan(1,length(CapRange));
SquaredSum_r0 = nan(1,length(CapRange));
Capped_intensity_profile_r0 = nan(41,1,length(CapRange));

for i=1:length(CapRange)
    % Cap the MS2 intensity with the value
    CapIntensity_r0(i) = CapRange(i);
    TotalmRNA_averaged_r0_Capped = TotalmRNA_averaged_r0;
    TotalmRNA_averaged_r0_Capped(TotalmRNA_averaged_r0 >= CapIntensity_r0(i)) = CapIntensity_r0(i);
    
    % Save the capped intensity profiles for a sanity check
     Capped_intensity_profile_r0(:,:,i) = TotalmRNA_averaged_r0_Capped;
     
    % Get the least square between in situ and MS2
    r0_insitu_norm = Averaged_Intensity_r0./max(Averaged_Intensity_r0);
    r0_MS2_norm = TotalmRNA_averaged_r0_Capped./ CapIntensity_r0(i);
    r0_MS2_norm(isnan(r0_MS2_norm)) = 0;
    
    % Get the squared sum for 20-60% of the embryo
    SquaredSum_r0(i) = sum((r0_insitu_norm(21:2.5:61) - r0_MS2_norm(9:25)').^2);
end

SquaredSum_r0(isnan(SquaredSum_r0)) = inf;
minIndex_r0 = find(SquaredSum_r0 == min(SquaredSum_r0));
cutOff_r0 = CapIntensity_r0(minIndex_r0);

% Sanity check
% hold on
% errorbar(AP, Averaged_Intensity_r0./max(Averaged_Intensity_r0), SEM_Intensity_r0./max(Averaged_Intensity_r0))
% for i=1:length(CapRange)
%     plot(0:0.025:1, Capped_intensity_profile(:,:,i)./max(Capped_intensity_profile(:,:,i)))
%     pause
% end

% r1
CapRange = 0:100:10000;

CapIntensity_r1 = nan(1,length(CapRange));
SquaredSum_r1 = nan(1,length(CapRange));
Capped_intensity_profile_r1 = nan(41,1,length(CapRange));

for i=1:length(CapRange)
    % Cap the MS2 intensity with the value
    CapIntensity_r1(i) = CapRange(i);
    TotalmRNA_averaged_r1_Capped = TotalmRNA_averaged_r1;
    TotalmRNA_averaged_r1_Capped(TotalmRNA_averaged_r1 >= CapIntensity_r1(i)) = CapIntensity_r1(i);
    
    % Save the capped intensity profiles for a sanity check
     Capped_intensity_profile_r1(:,:,i) = TotalmRNA_averaged_r1_Capped;
     
    % Get the least square between in situ and MS2
    r1_insitu_norm = Averaged_Intensity_r1./max(Averaged_Intensity_r1);
    r1_MS2_norm = TotalmRNA_averaged_r1_Capped./ CapIntensity_r1(i);
    r1_MS2_norm(isnan(r1_MS2_norm)) = 0;
    
    % Get the squared sum for 20-60% of the embryo
    SquaredSum_r1(i) = sum((r1_insitu_norm(21:2.5:61) - r1_MS2_norm(9:25)').^2);
end

SquaredSum_r1(isnan(SquaredSum_r1)) = inf;
minIndex_r1 = find(SquaredSum_r1 == min(SquaredSum_r1));
cutOff_r1 = CapIntensity_r1(minIndex_r1);

% r2
CapRange = 0:100:10000;

CapIntensity_r2 = nan(1,length(CapRange));
SquaredSum_r2 = nan(1,length(CapRange));
Capped_intensity_profile_r2 = nan(41,1,length(CapRange));

for i=1:length(CapRange)
    % Cap the MS2 intensity with the value
    CapIntensity_r2(i) = CapRange(i);
    TotalmRNA_averaged_r2_Capped = TotalmRNA_averaged_r2;
    TotalmRNA_averaged_r2_Capped(TotalmRNA_averaged_r2 >= CapIntensity_r2(i)) = CapIntensity_r2(i);
    
    % Save the capped intensity profiles for a sanity check
     Capped_intensity_profile_r2(:,:,i) = TotalmRNA_averaged_r2_Capped;
     
    % Get the least square between in situ and MS2
    r2_insitu_norm = Averaged_Intensity_r2./max(Averaged_Intensity_r2);
    r2_MS2_norm = TotalmRNA_averaged_r2_Capped./ CapIntensity_r2(i);
    r2_MS2_norm(isnan(r2_MS2_norm)) = 0;
    
    % Get the squared sum for 20-60% of the embryo
    SquaredSum_r2(i) = sum((r2_insitu_norm(21:2.5:61) - r2_MS2_norm(9:25)').^2);
end

SquaredSum_r2(isnan(SquaredSum_r2)) = inf;
minIndex_r2 = find(SquaredSum_r2 == min(SquaredSum_r2));
cutOff_r2 = CapIntensity_r2(minIndex_r2);

% r3
CapRange = 0:100:10000;

CapIntensity_r3 = nan(1,length(CapRange));
SquaredSum_r3 = nan(1,length(CapRange));
Capped_intensity_profile_r3 = nan(41,1,length(CapRange));

for i=1:length(CapRange)
    % Cap the MS2 intensity with the value
    CapIntensity_r3(i) = CapRange(i);
    TotalmRNA_averaged_r3_Capped = TotalmRNA_averaged_r3;
    TotalmRNA_averaged_r3_Capped(TotalmRNA_averaged_r3 >= CapIntensity_r3(i)) = CapIntensity_r3(i);
    
    % Save the capped intensity profiles for a sanity check
     Capped_intensity_profile_r3(:,:,i) = TotalmRNA_averaged_r3_Capped;
     
    % Get the least square between in situ and MS2
    r3_insitu_norm = Averaged_Intensity_r3./max(Averaged_Intensity_r3);
    r3_MS2_norm = TotalmRNA_averaged_r3_Capped./ CapIntensity_r3(i);
    r3_MS2_norm(isnan(r3_MS2_norm)) = 0;
    
    % Get the squared sum for 20-60% of the embryo
    SquaredSum_r3(i) = sum((r3_insitu_norm(21:2.5:61) - r3_MS2_norm(9:25)').^2);
end

SquaredSum_r3(isnan(SquaredSum_r3)) = inf;
minIndex_r3 = find(SquaredSum_r3 == min(SquaredSum_r3));
cutOff_r3 = CapIntensity_r3(minIndex_r3);

%% step2. Cap the amplitude of TotalmRNA using the CutOff value.
% First, define the values to be modified.

% Re-calculate the re-normalized r0,1,2,3 for MS2
TotalmRNA_averaged_r0_Capped(TotalmRNA_averaged_r0 >= cutOff_r0) = cutOff_r0;
TotalmRNA_averaged_r1_Capped(TotalmRNA_averaged_r1 >= cutOff_r1) = cutOff_r1;
TotalmRNA_averaged_r2_Capped(TotalmRNA_averaged_r2 >= cutOff_r2) = cutOff_r2;
TotalmRNA_averaged_r3_Capped(TotalmRNA_averaged_r3 >= cutOff_r3) = cutOff_r3;

TotalmRNA_SEM_r0_Capped = TotalmRNA_SEM_r0;
TotalmRNA_SEM_r1_Capped = TotalmRNA_SEM_r1;
TotalmRNA_SEM_r2_Capped = TotalmRNA_SEM_r2;
TotalmRNA_SEM_r3_Capped = TotalmRNA_SEM_r3;

TotalmRNA_SEM_r0_Capped(TotalmRNA_averaged_r0_Capped >= cutOff_r0) = 0;
TotalmRNA_SEM_r1_Capped(TotalmRNA_averaged_r1_Capped >= cutOff_r1) = 0;
TotalmRNA_SEM_r2_Capped(TotalmRNA_averaged_r2_Capped >= cutOff_r2) = 0;
TotalmRNA_SEM_r3_Capped(TotalmRNA_averaged_r3_Capped >= cutOff_r3) = 0;

%% step3. Plotting the re-normalized curves together

AP = 0.005:0.01:0.995;
APaxis = 0:0.025:1;

% r0
figure_compare_r0_CutOff = figure(1)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r0./max(Averaged_Intensity_r0), SEM_Intensity_r0./max(Averaged_Intensity_r0))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r0_Capped./cutOff_r0, TotalmRNA_SEM_r0_Capped./cutOff_r0)

title('Accumulated mRNA over AP - r0')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r1
figure_compare_r1_CutOff = figure(2)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r1./max(Averaged_Intensity_r1), SEM_Intensity_r1./max(Averaged_Intensity_r1))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r1_Capped./cutOff_r1, TotalmRNA_SEM_r1_Capped./cutOff_r1)
title('Accumulated mRNA over AP - r1')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r2
figure_compare_r2_CutOff = figure(3)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r2./max(Averaged_Intensity_r2), SEM_Intensity_r2./max(Averaged_Intensity_r2))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r2_Capped./cutOff_r2, TotalmRNA_SEM_r2_Capped./cutOff_r2)
title('Accumulated mRNA over AP - r2')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

% r3
figure_compare_r3_CutOff = figure(4)
hold on
% In Situ
errorbar(AP, Averaged_Intensity_r3./max(Averaged_Intensity_r3), SEM_Intensity_r3./max(Averaged_Intensity_r3))
% MS2-MCP
errorbar(APaxis, TotalmRNA_averaged_r3_Capped./cutOff_r3, TotalmRNA_SEM_r3_Capped./cutOff_r3)
title('Accumulated mRNA over AP - r3')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (Normalized)')
legend('In Situ','MS2')
StandardFigure(gcf,gca)

%% Plot CutOff values
cutOff_figure = figure;
plot([0 1 2 3], [cutOff_r0, cutOff_r1, cutOff_r2, cutOff_r3],'o')
xlim([-1 4])
ylim([0 10000])
title('Saturation values - optimized')
xlabel('Constructs')
ylabel('Saturation value')

StandardFigure(cutOff_figure, cutOff_figure.CurrentAxes)

%% Save plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InSitu_MS2_compare_plots\CutOff_Normalzied_profile_comparison_OptimizedCutOffs';

saveas(figure_compare_r0_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r0.tif'])
saveas(figure_compare_r0_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r0.pdf'])

saveas(figure_compare_r1_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r1.tif'])
saveas(figure_compare_r1_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r1.pdf'])

saveas(figure_compare_r2_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r2.tif'])
saveas(figure_compare_r2_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r2.pdf'])

saveas(figure_compare_r3_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r3.tif'])
saveas(figure_compare_r3_CutOff, [FigPath,filesep,'Comparison_TotlamRNA_InSitu_MS2_r3.pdf'])

saveas(cutOff_figure, [FigPath,filesep,'Comparison_SaturationPoints.tif'])
saveas(cutOff_figure, [FigPath,filesep,'Comparison_SaturationPoints.pdf'])


%% %%%%%%%%%%%%%%%%% Part2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Statistics for the boundary features
for i=1:length(ProcessedInSituData)
    BoundaryPosition(i) = ProcessedInSituData(i).BoundaryPosition;
    Slope(i) = ProcessedInSituData(i).Slope;
end

BoundaryPosition_mean = nanmean(BoundaryPosition);
BoundaryPosition_std = nanstd(BoundaryPosition);

%% r0 : assign the values to save
BoundaryPosition_mean_r0 = BoundaryPosition_mean;
BoundaryPosition_std_r0 = BoundaryPosition_std;

ProcessedInSituData_r0(1).BoundaryPosition_mean_r0 = BoundaryPosition_mean_r0;
ProcessedInSituData_r0(1).BoundaryPosition_std_r0 = BoundaryPosition_std_r0;
%% For r1, I chose the embryos (images) that are in NC14
% Image Name = Image[3,4,6,8,9,10,12,14]
% Index of those Images in the ProcessedInSituData
ProcessedInSituData_r1 = ProcessedInSituData;

Index_r1_NC14 = [2, 4, 6, 8, 9, 11, 13, 14];

% Initialize the Boundary Position/Slope 
BoundaryPosition_r1 = nan(1,length(ProcessedInSituData_r1));
Slope_r1 = nan(1,length(ProcessedInSituData_r1));

for i=Index_r1_NC14
    BoundaryPosition_r1(i) = ProcessedInSituData_r1(i).BoundaryPosition;
    Slope_r1(i) = ProcessedInSituData_r1(i).Slope;
end

BoundaryPosition_mean_r1 = nanmean(BoundaryPosition_r1);
BoundaryPosition_std_r1 = nanstd(BoundaryPosition_r1);

ProcessedInSituData_r1(1).Index_r1_NC14 = Index_r1_NC14;
ProcessedInSituData_r1(1).BoundaryPosition_mean_r1 = BoundaryPosition_mean_r1;
ProcessedInSituData_r1(1).BoundaryPosition_std_r1 = BoundaryPosition_std_r1;
%% For r2, I chose the embryos (images) that are in NC14
% Image Name = Image[15,21,23,25,28,29,31,33,34,35,38,39,42,43,44,45,46,51,52,53]
% Index of those Images in the ProcessedInSituData
ProcessedInSituData_r2 = ProcessedInSituData;

Index_r2_NC14 = [1,7,9,11,14,15,17,19,20,21,24,25,28,29,30,31,32,37,38,39];

% Initialize the Boundary Position/Slope 
BoundaryPosition_r2 = nan(1,length(ProcessedInSituData_r2));
Slope_r2 = nan(1,length(ProcessedInSituData_r2));

for i=Index_r2_NC14
    BoundaryPosition_r2(i) = ProcessedInSituData_r2(i).BoundaryPosition;
    Slope_r2(i) = ProcessedInSituData_r2(i).Slope;
end

BoundaryPosition_mean_r2 = nanmean(BoundaryPosition_r2);
BoundaryPosition_std_r2 = nanstd(BoundaryPosition_r2);

ProcessedInSituData_r2(1).Index_r2_NC14 = Index_r2_NC14;
ProcessedInSituData_r2(1).BoundaryPosition_mean_r2 = BoundaryPosition_mean_r2;
ProcessedInSituData_r2(1).BoundaryPosition_std_r2 = BoundaryPosition_std_r2;

end