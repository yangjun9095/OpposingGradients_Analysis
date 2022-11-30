% Sequence of pipeline scripts to run for the image analysis of OpposingGradients datasets
% Last updated : 7/22/2020 by Yang Joon Kim
%

%% Export the LIF file
% Note that now we don't have to export the .lif files to .tif files with
% LASX.
Prefix = ExportDataForLivemRNA
% Click the folder containing the raw movie (.lif file, .lif file does not
% have to be exported.)
% Make sure to have a good nuclear channel. Usually, the default option for
% z-projection works very well, so no need to change.
% The cool thing is that you can define mitoses in this step in the GUI.

%% filtering movie
filterMovie(Prefix,'highPrecision','customFilter','Difference_of_Gaussian_3D','nWorkers',1)

%% spot segmentation
% Thresh = ? % after ImageJ looking through z-stacks and frames.
segmentSpots(Prefix, Thresh,'Shadows',1,'keepProcessedData')%'fit3D','displayFigures')

%% TrackNuclei
TrackNuclei(Prefix)

%% CheckNucleiSegmentation
CheckNucleiSegmentation(Prefix)
%% TrackmRNADynamics
TrackmRNADynamics(Prefix)

%% CheckParticleTracking 
CheckParticleTracking(Prefix)
% Let's make sure to check this before we launch into the next steps.
% Reminder : This step is checking 
% 1) whether we have a good spot segmentation.
% 2) whether we have a good particle tracking.
% language : We use "Spots" for spots segmented for each frame, then we use
% "Particles" for the temporally tracked spots.
% Let's check the two things above for most of the NC14 frames. GO through
% some single particle time trace by clicking "p" then moving using "k" is
% also a good idea. If the particle tracking was done well, you should be
% able to see a nice trace that looks like averaged time traces.

%% FindAPAxisFullEmbryo
% AP axis registration
FindAPAxisFullEmbryo(Prefix)

%% AddParticlePosition
AddParticlePosition('ManualAlignment','SelectChannel')
% Make sure to select the Histone channel as the marker.

%% CompileParticle
CompileParticles(Prefix,'skipAll','ApproveAll','MinParticles',2,'MinTime',2)

% Note that if we have a perfect particle curation, then we'd not need
% 'ApproveAll' option in here. 