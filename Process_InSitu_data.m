function [Intensity_raw,...
            nBins, Intensity_mean,  Intensity_std, ...
            window_smooth, Intensity_mean_smooth,...
            background_intensity, Intensity_mean_smooth_BGsubtracted] = Process_InSitu_data(InSituEmbryo_raw, EmbryoName,filePath)
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
Path = [filePath,filesep,'APdefined']
saveas(gcf,[Path,filesep,EmbryoName])
clf
%% Step2. Define the box window 
% Let's try to use the improfile
Intensity_raw = improfile(InSituEmbryo,X,Y)

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
background_intensity = nanmean(Intensity_mean_smooth(71:100));

Intensity_mean_smooth_BGsubtracted = Intensity_mean_smooth - background_intensity;

% For the negative intensity, I'll assume that they are actually zeros, and
% the negative value comes from the measurement noise, etc.
Intensity_mean_smooth_BGsubtracted(Intensity_mean_smooth_BGsubtracted<0) =0;

end