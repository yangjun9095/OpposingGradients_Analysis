function [Intensity_raw,...
            nBins, Intensity_mean,  Intensity_std, ...
            window_smooth, Intensity_mean_smooth,...
            background_intensity, Intensity_mean_smooth_BGsubtracted, X, Y] = Process_InSitu_data(InSituEmbryo_raw, EmbryoName,filePath)
%% This function grabs the raw .jpg file of In Situ embryo image, then process it 
% 1) We can define the AP axis by clicking two points
% 2) We can bin it with 100 AP bins, also smoothen it (this smoothing
% window could be fed into as an option.)
% 3) Background subtraction
%% Step0. Convert RGB file to grayscale (8-bit), then invert the contrast

InSituEmbryo_8bit = rgb2gray(InSituEmbryo_raw);

InSituEmbryo = imcomplement(InSituEmbryo_8bit);



%% Step1. Define the AP axis
% This is done by clicking two points, using ginput

% % For now, let's just assume that the embryos are alingned perfectly
% % horizontal, thus I can just get the anterior point, and draw a horizontal
% % line.

% This has to be re-done later by setting up the AP line, then finding the pixels
% based on that.
hold on
imshow(InSituEmbryo)

% Define the AP axis by clicking the head and tail points.
[X, Y] = ginput(2);

head_x = X(1);
head_y = Y(1);

tail_x = X(2);
tail_y = Y(2);

hold on
plot(head_x,head_y,'g.','MarkerSize',20);
plot(tail_x,tail_y,'r.','MarkerSize',20);
plot(X,Y,'k')

mkdir([filePath,filesep,'APdefined'])
mkdir([filePath,filesep,'ProcessedResults'])
Path = [filePath,filesep,'APdefined'];
saveas(gcf,[Path,filesep,EmbryoName])
close gcf
%% Step2. Define the box window 
% Let's try to use the improfile
Width = 30; % 30 pixels used in Chen, 2012 for the width
[CX,CY,C_sum,C,xi,yi,im_pi2] = improfile_integrated(InSituEmbryo,Width,X,Y);
Intensity_raw = C_sum/Width;
% Intensity_raw_single = improfile(InSituEmbryo,X,Y);

% In case the image is fipped (left to right, then flip the intensity
% vector)
if Intensity_raw(200) < Intensity_raw(end-200) 
    Intensity_raw = fliplr(Intensity_raw)
    disp('Flip')
end

%% Step3. Averaging over some window
% Here, I need a good way of binning, for example, like in 100 bins as
% Jeehae did.
nBins = 100;
deltaEL = 1/nBins;
EL = 0:deltaEL:1;

lengthAP = length(Intensity_raw);

binSize = lengthAP/nBins;

binFilter = floor((1:100) * binSize);
binFilter = [1 binFilter];

for i=1:nBins
    Intensity_mean(i) = nanmean(Intensity_raw(binFilter(i:i+1)));
    Intensity_std(i) = nanstd(Intensity_raw(binFilter(i:i+1)));
end

%% Optional : Smoothening the curve
window_smooth = 5;
Intensity_mean_smooth = movmean(Intensity_mean,window_smooth);

%% Step5. (Optional?) Background subtraction
% The issue right now is that the background is underestimated when we only
% think about the minimum value across AP.
% Here, I'll make an assumption that there's no expression after 70% of the
% embryo length (AP axis). 
background_intensity = nanmean(Intensity_mean_smooth(71:90));

Intensity_mean_smooth_BGsubtracted = Intensity_mean_smooth - background_intensity;

% For the negative intensity, I'll assume that they are actually zeros, and
% the negative value comes from the measurement noise, etc.
Intensity_mean_smooth_BGsubtracted(Intensity_mean_smooth_BGsubtracted<0) =0;

%% Plot the processed In Situ intensity profile 
% %(optional, this should be done in the upstream script, which is Supple_02_compare_InSitu_MS2MCP.m
% 
% % Define the Figure path
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\InSitu_MS2_compare_plots';
% 
% hold on
% % Raw intensity
% plot((1:length(Intensity_raw))/length(Intensity_raw) ,Intensity_raw)
% % Averaged pixel intensity, std (over AP)- binned
errorbar(0.005:0.01:0.995,Intensity_mean,Intensity_std)
close gcf
% % Smoothened over window_smooth
% plot(0.005:0.01:0.995,Intensity_mean_smooth)
% 
% % plot(BoundaryPosition, Intensity_mean_smooth(AP == BoundaryPosition),'o','MarkerSize',10)
% % plot(BoundaryPosition_inflection,Intensity_mean_smooth(AP == BoundaryPosition_inflection),'o','MarkerSize',10)
% 
% % Bar line for the Background
% % plot(AP, ones(size(AP))*Background,'k')
% 
% % tangential line
% % plot(AP,Slope*(AP-BoundaryPosition) + Intensity_mean_smooth(AP == BoundaryPosition))
% 
% ylim([0 max(Intensity_mean_smooth)+50])
% title('Intensity from In Situ embryo')
% xlabel('AP (EL)')
% ylabel('Intensity (AU)')
% legend('Raw','Binned','Smoothened')
% StandardFigure(gcf,gca)
% 
% % Save the result figure
% saveas(gcf,[FigPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'.tif'])
% saveas(gcf,[FigPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'.pdf'])
% 
% %% Optional : Background-subtracted profile
% plot(0.005:0.01:0.995,Intensity_mean_smooth_BGsubtracted)
% title('Intensity from In Situ embryo-BG subtracted')
% xlabel('AP (EL)')
% ylabel('Intensity (AU)')
% legend('In Situ intensity')
% StandardFigure(gcf,gca)
% 
% % Save the result figure
% saveas(gcf,[FigPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'_BGsubtracted.tif'])
% saveas(gcf,[FigPath,filesep,'ProcessedResults',filesep,EmbryoName(1:end-4),'_BGsubtracted.pdf'])

end