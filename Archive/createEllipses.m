function createEllipses(Prefix)
%% Description
% This is a function to create a fake Ellipses.mat structure
% such that we can add nuclei at certain frames by manually, then copy,
% paste for the rest of the frames.
% For convenience, I'll add two nuclei per frame so that the code can run
% without a bug

% Prefix
% PreProc, read His images
% Create 

%% Define the file paths
% [~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders;
% Define the Dropbox folder path as in OpposingGradientsResults
DropboxFolder = 'S:\YangJoon\Dropbox\OpposingGradient';

% Import the FrameInfo to determine the length of the movie
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'])
nFrames = length(FrameInfo);
%% Structure of an Ellipses.mat file
 load(['S:\YangJoon\Dropbox\OpposingGradient\2019-10-13-hbP2-r2_far_MS2V5-lacZ-4',filesep,'2019-10-13-hbP2-r2_far_MS2V5-lacZ-4_lin.mat'])
% Cell of (number of time frames x 1)
% For each cell : (number of nucleus x 9) double
% For each nucleus, there are xPos, yPos, xDiameter, yDiameter, etc.

% Create an empty cell whose size is the Ellipses.mat

Ellipses = cell(nFrames,1);
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses');
end