function FitSingleParticleLinearRise(varargin)

%Last updated on: 3/2/2020 by Jonathan Liu

%
%
%Fits linear rise to steady state for single particle MS2 signals, to
%extract a t_on at the single cell level. This script only looks at nc13
%for now.

%Variable inputs (Name/Value pairs):
%   'Prefix', Prefix: Input a Prefix string. If not chosen, user
%                           selects Prefix in a dialog box.
%   'NCWindow', [nc13start, nc13end]: Start and end
%                           times of fitting for cycles 13 and 14.
%                           Pass this in as a 1x2 array of doubles.
%                           Read as [time after nc13 start, time after nc14 start].
%                           Default values are [1, -5].
%   'KeepPool', []: Keeps the parallel pool after running (no 2nd argument
%                   required)

%% Variable inputs

%Default settings
LoadPrefix = true; %User selection of Prefix by default.
ncwindow = [1, -4]; %Time of MCMC fitting for nc13 and nc14
KeepPool = false; %By default, shut down parallel pool after running code

%User specified settings
for i=1:length(varargin)
    if strcmpi(varargin{i},'Prefix')
        Prefix = varargin{i+1};
        LoadPrefix = false;
    end
    if strcmpi(varargin{i},'NCWindow')
        ncwindow = varargin{i+1};
    end
    if strcmpi(varargin{i},'KeepPool')
        KeepPool = true;
    end
end

%Extract nuclear cycle fit windows.
firstnc13time = ncwindow(1,1);
lastnc13time = ncwindow(1,2);

%% Load data
%Get the default folders for Prefix loading
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;

if LoadPrefix
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
end

%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);
        
%Load the compiled particles and the division information                                    
data = load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat']);

% Extract the fields from the cell structure (This is for fields like
% CompiledParticles that are saved inside {}.
channel = 1; %We're working with one-color data for now.

if iscell(data.CompiledParticles)
    CompiledParticles = data.CompiledParticles{channel};
else
    CompiledParticles = data.CompiledParticles; %Legacy data format compability
end

%% Set some initial definitions

%Definitions
nc13 = data.nc13; %nc13 frame index
nc14 = data.nc14; %nc14 frame index
t = data.ElapsedTime - data.ElapsedTime(nc13); %time variable for nc13
N_particles = length(CompiledParticles); %total number of particles

%% Loop through single particles for fit

clear FitResults
clear ResultsPlot

FitResults = [];
ResultsPlot = [];

w = waitbar(0,'Fitting nuclei...');
for nuc = 1:N_particles
waitbar(nuc/N_particles,w,['Fitting nucleus ',num2str(nuc),' of ',num2str(N_particles)]);

%MS2 fluorescence data
Frame_raw = CompiledParticles(nuc).Frame;
MS2_raw = CompiledParticles(nuc).Fluo;

%Remove datapoints before and after the window of fitting
FramesToKeep = and(t(Frame_raw) > t(nc13) + firstnc13time, t(Frame_raw) < t(nc14) + lastnc13time);
Frame = Frame_raw(FramesToKeep);
MS2 = MS2_raw(FramesToKeep);
timepoints = t(Frame); %Timepoints of MS2 signal

%Skip this nucleus if there are fewer than 10 datapoints in the nuclear
%cycle
Threshold = 10;
if length(Frame) < Threshold
    disp(['skipping particle ',num2str(nuc),' of ',num2str(N_particles), ' (no data)']);
    continue
end
disp(['Fitting particle ',num2str(nuc),' of ',num2str(N_particles)]);
%% Run single particle fit
%Make interpolated time vector for fitting script
t_fold = 10; %For now, interpolate time with 10x resolution to get a better resolution on t_on.
dt = t(2) - t(1); %Experimental time resolution
t_interp = t(nc13):(dt/t_fold):t(nc14); %Fit time vector

%Get fit parameters
[basal, t_on, t_steady, rate, resnorm] = FitLinearRise(MS2,timepoints,t_interp);

%Plotting variables
t_fit = t_interp;
LinearFit = MakeLinearRise(basal,t_on,t_steady,rate,t_interp);
t_plot = timepoints;
MS2_plot = MS2;

%% Save data
FitResults(end+1).basal = basal;
FitResults(end).t_on = t_on;
FitResults(end).t_steady = t_steady;
FitResults(end).rate = rate;
FitResults(end).resnorm = resnorm;
FitResults(end).MeanAP = CompiledParticles(nuc).MeanAP;
FitResults(end).ApprovedFits = 0; %Set approval status to uncurated (0).

ResultsPlot(end+1).t_plot = t_plot;
ResultsPlot(end).MS2_plot = MS2_plot;
ResultsPlot(end).t_fit = t_fit;
ResultsPlot(end).LinearFit = LinearFit;
ResultsPlot(end).firstnc13time = firstnc13time;
ResultsPlot(end).lastnc13time = lastnc13time;
end

%% Save data into .mat structure
loc = [DropboxFolder,filesep,Prefix,'\'];
filename = 'LinearRiseFits';
save([loc,filename,'.mat'],'FitResults','ResultsPlot','Prefix');

close(w);
disp(['Linear fit analysis complete. Information stored in: ',loc]);

%Close the parallel pool if desired
if ~KeepPool
    poolobj = gcp('nocreate');
    delete(poolobj);
end
end

function [basal, t_on, t_steady, rate, resnorm] = FitLinearRise(MS2,timepoints,t_interp)
%Fitting script to fit linear rise to MS2 signal

%Anonymous fitting function

%Interpolate signal to data
fun = @(x)interp1(t_interp,MakeLinearRise(x(1),x(2),x(3),x(4),t_interp),timepoints) - MS2;

%Fit the data
basal0 = 1;
t_on0 = 5;
t_steady0 = 10;
rate0 = 10;
x0 = [basal0, t_on0, t_steady0, rate0]; %Initial condition
lb = [0, 0, 0, 0]; %lower bound
ub = [1000, 20, 20, 1000]; %upper bound
options = optimoptions(@lsqnonlin,'Display','off');

[x,resnorm] = lsqnonlin(fun,x0,lb,ub,options); %Fit

%Extract fit parameteres
basal = x(1);
t_on = x(2);
t_steady = x(3);
rate = x(4);
resnorm = resnorm ./ length(MS2); %Normalize residual norm by number of datapoints.

end

function signal = MakeLinearRise(basal,t_on,t_steady,rate,t_interp)
%Function to generate a linear rise in signal to steady state
signal = zeros(size(t_interp)); %Output signal

%If t_steady is less than t_on, set them equal to enforce a positive signal
if t_steady <= t_on
    t_steady = t_on;
end

%Generate signal
for i = 1:length(t_interp)
    if t_interp(i) > t_on && t_interp(i) < t_steady
        signal(i) = signal(i) + rate*(t_interp(i) - t_on);
    elseif t_interp(i) > t_on && t_interp(i) >= t_steady
        signal(i) = signal(i) + rate*(t_steady - t_on);
    end
end

%Add basal fluorescence
signal(signal < basal) = basal;
end