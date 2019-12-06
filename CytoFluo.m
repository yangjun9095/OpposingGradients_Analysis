function [CytoFluoDensity] = CytoFluo(varargin)
% This function calculates the average fluorescence intensity of the whole cytoplasm for
% each frame and z slice. 
% It does so by first creating a binary mask using the nuclear segmentation information stored in Ellipses.mat
% it then applies this mask to the PreProcessed Bcd files

% Input: A movie prefix

% Output: a [frame x Z] array containing the average pixel intensity of the
% cytoplasm for each z slice in each frame.

% Authors: Simon Alamos and Jordan Xiao (modified by YangJoon Kim)
% Contact: simon.alamos@berkeley.edu

% November 2017
%% Get Folder info, set up stuff, etc

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Prefix=varargin{1};

FilePrefix=[Prefix,'_'];

%Now get the actual Dropbox folder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

% This is to get the name of the channels from MovieDatabase
[XLSNum,XLSTxt,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);

%Load all the information

%load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
load([DropboxFolder,filesep,Prefix,filesep,[Prefix '_lin.mat']])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])


%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
end

% %get frames of each mitosis
[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
% nc10Column=find(strcmp(XLSRaw(1,:),'nc10'));
% nc11Column=find(strcmp(XLSRaw(1,:),'nc11'));
% nc12Column= strcmp(XLSRaw(1,:),'nc12');
% nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
% nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

%% Create a 3D mask using ellipses

% Our nucleus segmentation algorithm  DOES NOT track nuclei in Z. In
% consequence, a single segmentation is applied to all Z slices in any
% given frame. This means we don't need 4 dimensions, but just 3: X,Y,T

% create a black image of the same X,Y size as the data

%get dimensions
TotalFrames = length(FrameInfo);
Rows = FrameInfo(1).LinesPerFrame;
Columns = FrameInfo(1).PixelsPerLine;

%create black image
BlackImage = zeros(Rows,Columns);

%create empty XYT array to store masks
MovieMask = ones(Rows,Columns,TotalFrames);


figure(1)
h = waitbar(0,'Please wait, creating a XYT 3D mask for the cytoplasm');
for frame = 1:TotalFrames
    
    Mask = BlackImage; %we will add this frame's ellipses to this black image
    
    FrameEllipses = Ellipses{frame}; %get the ellipses belonging to this frame
    CentersX = FrameEllipses(:,1); %get the centers of this frame's ellipses
    CentersY = FrameEllipses(:,2); %get the sizes of this frame's ellipses
    EllipsesRadius = FrameEllipses(1,3); %get the radius of the ellipses in this frame
    
    %now we will add the ellipses one by one as white circles
    for ellipse = 1:size(FrameEllipses,1)
        
        % get parameters needed to create a disk
        SingleEllipseMask = zeros(Rows,Columns);
        EllipseXCenter = CentersX(ellipse);
        EllipseYCenter = CentersY(ellipse);
        Radius = 1.45*(ceil(EllipsesRadius)); %radius of the ellipse. Increase by 45% to make sure we get the whole nucleus
        
        
        % create a white disk where the ellipse is
        Circle=BlackImage;
        Circle=MidpointCircle(Circle,Radius,EllipseYCenter,...
        EllipseXCenter,1); 
    
        % add disk (ellipse) to the mask
        Mask = Mask+Circle; 

    end
    
    Mask = Mask==0; % invert the mask
    MovieMask(:,:,frame) = Mask;  %add mask to the XYT movie array of masks
    imshow(Mask,[])
    %pause(0.01)
    waitbar(frame / TotalFrames);


end
close(h)
%% Apply mask to Bcd / Dl channel movie

% First, we will create a 4D movie using the Bcd / Dorsal channel so that we can
% apply the mask to it by doing element-wise matrix multiplication.

%get the size of the Z dimension
Slices = FrameInfo(1).NumberSlices + 2;

% Figure out which channel has Bcd or Dorsal
if contains(lower(Channel1{1}),'dorsal') || contains(lower(Channel1{1}),'dl')||contains(lower(Channel1{1}),'bcd')
    InputChannel = 'ch01';
elseif contains(lower(Channel2{1}),'dorsal') || contains(lower(Channel2{1}),'dl')||contains(lower(Channel2{1}),'bcd')
    InputChannel = 'ch02';
else
    display ('no channel with Bcd or Dl in its name found in this dataset')
end

InputChannel = 'ch01';

InputChannelFiles = dir([PreProcPath,filesep,Prefix,filesep,'*' InputChannel '.tif']); % save info about files in a struct
Frames = length(InputChannelFiles)/Slices; % get the number of frames

%%  We need add to the struct the Z and T info of each file, which is stored in their names

for i= 1:length(InputChannelFiles)
    
    ImageName = InputChannelFiles(i).name; %get the name string
    SplitName = strsplit(ImageName,'_'); %split the name in strings between '_' characters
    
    % Z information
    ZSlice1 = SplitName{end-1}; %select only the part of the name referring to the Z info
    ZSlice2 = str2num(ZSlice1(2:end)); %get an actual number corresponding to the z slice
 
    % Time information
    Time1 = SplitName{end-2}; %select only the part of the name referring to the T info
    Time2 = str2num(Time1); %get an actual number corresponding to the frame
    
    % Add info to the struct
    InputChannelFiles(i).Z = ZSlice2;
    InputChannelFiles(i).Time = Time2;
    
end

%% Now, loop over time and z to apply the cytoplasm mask to each image
% We will store the cytoplasm fluorescence in a T x Z 2D matrix called
% 'CytoFluo'
% For each frame we will calculate the mean cytoplasmic fluorescence of each Z slice.

CytoFluoDensity = zeros(Frames,Slices); %initialize the array where we'll store the data
MovieMask(MovieMask==0)=NaN; %convert 0s to NaNs because the raw data contains 0s

for T = 1:Frames
    
    Mask = MovieMask(:,:,T);
    
    
    for Z = 1:Slices
        
        for i = 1:length(InputChannelFiles)
            
            ImageName = InputChannelFiles(i).name;
            ZSlice = InputChannelFiles(i).Z;
            Time = InputChannelFiles(i).Time;
            
            if ZSlice == Z && Time == T
                
                NucleusFluoImage = imread(ImageName); % read the image
                % imultiply sometimes doesn't work, depending on the
                % matlab version because of the class mismatch
                NucleusFluoImage = double(NucleusFluoImage);
                CytoFluoImage = immultiply(NucleusFluoImage,Mask); % apply mask to image
                %imshow(CytoFluoImage,[])
                %title([num2str(i) num2str(Z)])
                %waitforbuttonpress
                
                % save the mean cytoplasmic fluorescence of this Z slice in
                % this frame
                CytoFluoDensity(T,Z) = nanmean(CytoFluoImage(:));              
                
            end
        end
    end
end

% save results in the DynamicsResults folder
save([DropboxFolder,filesep,Prefix,filesep,'MeanCytoFluo.mat'],'CytoFluoDensity')

%% Now let's actually caclulate the difference in fluorescence between nucleus and cytoplasm
% this part was adapted from Jordan's script

% We need to get the integration area (in pixels) that was used to
% calculate the nuclear fluorescence density so that the units match with the cyto
% fluo denisty. This corresponds to the variable 'IntegrationRadius', which
% is defined in the function 'TrackNuclei' as being 2 (in microns)

absIntegrationRadius = 2; %in microns
pixelSize = FrameInfo(1).PixelSize; %in microns/pixel

% divide the standard two micrometers (assumed actual nucleus size) by the resolution ([pixelsize] = micrometers per pixel)
IntegrationRadius = floor(absIntegrationRadius/pixelSize); %in pixels
% not sure what this is for but it's in TrackNuclei so...
if ~mod(IntegrationRadius,2)
    IntegrationRadius=IntegrationRadius+1;
end

% This is taken from ExtractNuclearFluorescence. We use this to
% calculate the integration area.
Circle = logical(zeros(3*IntegrationRadius,3*IntegrationRadius));
Circle = MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1);
integrationArea = sum(sum(Circle)); %number of pixels in the integration area

% Loop over schnitzcells to calculate the delta (Nucleus-Cyto) for each one
%for each schnitz, divide every fluorescence in schnitzcells by the integration area
% for s = 325:330%length(schnitzcells)
%     s
%     timeframes = schnitzcells(s).frames; % vector with all the timeframes for a given schnitz
%     fluos = schnitzcells(s).Fluo; % array with one row per timeframe, one column per z slice
%     fluodens = fluos/integrationArea; % average fluorescence intensity per pixel, to compare with CytoFluoDens
%     dFluo = zeros(size(fluos)); % vector to be added to final cell array
% 
%     for f = 1:length(timeframes) % t is the index of the frame within timeframes vector, NOT the real timeframe
%         frame = timeframes(f); % the real timeframe within the movie
%         for z = 2:size(fluodens,2)-1 % z index--ignore top and bottom artificial slices
%             nucFluo = fluodens(f,z);
%             cytoFluo = CytoFluoDensity(frame,z);
%             dFluo(f,z) = abs(nucFluo - cytoFluo); % absolute difference between cyto and nuclear fluo densities
%         end
%         % plot delta
%         plot(dFluo(f,:),'k')
%         hold on
%         % plot nucleus
%         plot(fluodens(f,:),'r')
%        % plot cytoplasm
%         plot(CytoFluoDensity(f,:),'b')
%         hold off
%         legend('Delta','Nucleus','Cytoplasm')
%         waitforbuttonpress  
%     end
% 
%     %add the delta fluo to the struct
%     schnitzcells(s).DeltaFluo = dFluo;
% 
%     end
end

    


