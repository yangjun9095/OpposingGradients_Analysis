function Supple_02_compare_InSitu_MS2MCP
% DESCRIPTION
% This script is for comparison of In Situ hybridization and MS2-MCP 
% in terms of accumulated cytoplasmic mRNA pattern. 
% As the first step, we will compare Steve Small lab's data for P2-r0,1,2,3, with our
% P2(r0,1,2,3)-MS2-MCP data.

%% Load the In Situ datasets
InSituPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Embryos for Yang Joon';
r0Path = [InSituPath,filesep,'hbP2\20140303 G49 hb-lacZ 1-1\NC14'];
r1Path = [InSituPath,filesep,'hbP2+1Run\20140304 cm4A hb-lacZ 1-1'];
r2Path = [InSituPath,filesep,'hbP2+2Run\20140304 cm2A hb-lacZ 1-1'];
r3Path = [InSituPath,filesep,'hbP2+2Run\20140304 cm2A hb-lacZ 1-1'];

%% Load specific embryos
% Use dir function to read out all names in the directory.
D = dir(r2Path);

%% Loop through all images (embryos)

% Define which files to import (Option)
% ProcessedInSituData_r0 = load([r0Path,filesep,'ProcessedInSituData_r0.mat']);
% ProcessedInSituData_r0 = ProcessedInSituData_r0.ProcessedInSituData_r0;

DataPath = r2Path;
%k=1; % Count the number of images

AP = 0.005:0.01:0.995;
    
for i=32:35%3:length(D)
    % import the raw image
    EmbryoName = D(i).name;
    InSituEmbryo_raw = imread([DataPath,filesep, EmbryoName]);
    % Get the raw, binned, smoothened, and backgroud-subtracted intensity
    % profile over the defined AP axis.
    [Intensity_raw,...
            nBins, Intensity_mean,  Intensity_std, ...
            window_smooth, Intensity_mean_smooth,...
            background_intensity, Intensity_mean_smooth_BGsubtracted] =...
                        Process_InSitu_data(InSituEmbryo_raw, EmbryoName, DataPath)
        
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
    
    %k=k+1;
end

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
    %saveas(gcf,[DataPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'.pdf'])
end

%% Get the statistics

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
%% Save the fields
% I want to save the processed data
% 1) the raw embryo image with AP axis defined (DONE)
% 2) Intensity (raw, binned, smoothened)
% 3) Boundary Position, slope
%ProcessedInSituData_r1 = ProcessedInSituData;

save([r0Path,filesep,'ProcessedInSituData_r0.mat'],'ProcessedInSituData_r0')
save([r1Path,filesep,'ProcessedInSituData_r1.mat'],'ProcessedInSituData_r1')
save([r2Path,filesep,'ProcessedInSituData_r2.mat'],'ProcessedInSituData_r2')
end