function FitMeanVectorWholeCycleNC14(varargin)

%Last updated on: 4/16/2020 by Jonathan Liu

%
%
%Fits whole nuclear cycle for mean MS2 signals, to
%extract various parameters binned by AP position. This script only looks at nc14.

%Variable inputs (Name/Value pairs):
%   'Prefix', Prefix: Input a Prefix string. If not chosen, user
%                           selects Prefix in a dialog box.
%   'NCWindow', [nc14start, nc14end]: Start and end
%                           times of fitting for cycle 14.
%                           Pass this in as a 1x2 array of doubles.
%                           Read as [time after nc14 start, time after nc14 start].
%                           Default values are [1, 25].
%   'LocalMovieDatabase', LocalDropboxFolderString: uses a local MovieDatabase
%   file contained in alternate DropboxFolder than default.
%   LocalDropboxFolderString is the string specifying the name of the
%   DropboxFolder as specified in ComputerFolders.

%% Variable inputs

%Default settings
LoadPrefix = true; %User selection of Prefix by default.
ncwindow = [1, 25]; %Time of fitting for nc13 and nc14
KeepPool = false; %By default, shut down parallel pool after running code

%User specified settings
UseLocalMovieDatabase = false; %By default, don't use a local MovieDatabase
for i=1:length(varargin)
    if strcmpi(varargin{i},'Prefix')
        Prefix = varargin{i+1};
        LoadPrefix = false;
    end
    if strcmpi(varargin{i},'NCWindow')
        ncwindow = varargin{i+1};
    end
    if strcmpi(varargin{i},'LocalMovieDatabase')
        UseLocalMovieDatabase = true;
        LocalDropboxFolderString = varargin{i+1};
    end
end

%Extract nuclear cycle fit windows.
firstnc14time = ncwindow(1,1);
lastnc14time = ncwindow(1,2);

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

% Extract the fields from the cell structure (This is for fields like MeanVectorAP
% that are saved inside {}.
channel = 1; %We're working with one-color data for now.

if iscell(data.MeanVectorAP)
    MeanVectorAP = data.MeanVectorAP{channel};
    SDVectorAP = data.SDVectorAP{channel};
    NParticlesAP = data.NParticlesAP{channel};
else
    MeanVectorAP = data.MeanVectorAP;
    SDVectorAP = data.SDVectorAP;
    NParticlesAP = data.NParticlesAP;
end

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DropboxFolder);

%% Set some initial definitions

%Definitions
IncompleteNC14Flag = false; %Flag if start of nc14 wasn't captured
nc14 = data.nc14; %nc14 frame index

%If start of nc14 wasn't captured, skip the fits
if nc14 ==0
    IncompleteNC14Flag = true;
    nc14 = 1;
    disp('Start of nc14 was not captured, skipping this dataset');
end
t = data.ElapsedTime - data.ElapsedTime(nc14); %time variable for nc14
bins = 0:APResolution:1; %AP Resolution
N_bins = length(bins); %Number of AP positions

%% Loop through AP positions for fit

clear FitResults
clear ResultsPlot

FitResults(N_bins) = struct;
ResultsPlot(N_bins) = struct;

%Initialize approval status and AP position
for AP = 1:N_bins
    FitResults(AP).ApprovedFits = 0;
    FitResults(AP).APPos = bins(AP);
end

w = waitbar(0,'Fitting nuclei...');
for AP = 1:N_bins
waitbar(AP/N_bins,w,['Fitting AP position ',num2str(AP),' of ',num2str(N_bins)]);

if IncompleteNC14Flag
    continue;
end

%Flag for automatic rejection for certain filter criteria
rejectApprovedFits = false;

%Nucleus index
APPos = bins(AP);

%MS2 fluorescence data
Frame_raw = 1:length(t); %Frame indices
MS2_raw = MeanVectorAP(:,AP)'; %Raw MS2 signal for this AP bin
SD_MS2_raw = SDVectorAP(:,AP)'; %Standard error in MS2 signal
NParticles_raw = NParticlesAP(:,AP)'; %Number of particles in this AP bin

%Remove datapoints before and after the window of fitting
FramesToKeep = and(t(Frame_raw) > t(nc14) + firstnc14time, t(Frame_raw) < t(nc14) + lastnc14time);
Frame = Frame_raw(FramesToKeep);
MS2 = MS2_raw(FramesToKeep);
SD_MS2 = SD_MS2_raw(FramesToKeep);
NParticles = NParticles_raw(FramesToKeep);

%Remove nan values
nanInd = isnan(MS2); %Nan indices
Frame(nanInd) = [];
MS2(nanInd) = [];
SD_MS2(nanInd) = [];
NParticles(nanInd) = [];
timepoints = t(Frame); %Timepoints of MS2 signal

%Skip this AP position if there are fewer than 10 datapoints in the nuclear
%cycle
Threshold = 10;
if sum(~isnan(MS2)) < Threshold
    disp(['skipping AP position ',num2str(AP),' of ',num2str(N_bins), ' (no data)']);
    continue
end
disp(['Fitting AP position ',num2str(AP),' of ',num2str(N_bins)]);
%% Run fit
%Make interpolated time vector for fitting script
t_fold = 10; %For now, interpolate time with 10x resolution to get a better resolution on t_on.
dt = t(2) - t(1); %Experimental time resolution
t_interp = t(nc14):(dt/t_fold):t(end); %Fit time vector

%Get fit parameters
[basal,t_on,t_dwell,t_active,tau,rate,residual,resnorm,CI] = FitWholeCycle(MS2,timepoints,t_interp);

%Plotting variables
t_fit = t_interp;
WholeCycleFit = MakeWholeCycle(basal,t_on,t_dwell,t_active,tau,rate,t_interp);
t_plot = timepoints;
MS2_plot = MS2;
SD_MS2_plot = SD_MS2;
NParticles_plot = NParticles;

%Filters to set approval status a priori
N_linear = sum(and(t_plot > t_on, t_plot < (t_on + t_dwell))); %Number of datapoints in linear regime
if N_linear < 3 %Require at least three datapoints to be in linear regime
    rejectApprovedFits = true;
end

%% Save data
FitResults(AP).basal = basal;
FitResults(AP).t_on = t_on;
FitResults(AP).t_dwell = t_dwell;
FitResults(AP).t_active = t_active + t_dwell; %Note that fitted t_active is measured from after the dwell time
FitResults(AP).tau = tau;
FitResults(AP).rate = rate;
FitResults(AP).residual = residual;
FitResults(AP).resnorm = resnorm;
FitResults(AP).CI = CI;
FitResults(AP).APPos = APPos;

%Automatically reject if not enough datapoints were in the linear regime
if rejectApprovedFits
    FitResults(AP).ApprovedFits = -1;
end

ResultsPlot(AP).t_plot = t_plot;
ResultsPlot(AP).MS2_plot = MS2_plot;
ResultsPlot(AP).SD_MS2_plot = SD_MS2_plot;
ResultsPlot(AP).NParticles_plot = NParticles_plot;
ResultsPlot(AP).t_fit = t_fit;
ResultsPlot(AP).WholeCycleFit = WholeCycleFit;
ResultsPlot(AP).firstnc14time = firstnc14time;
ResultsPlot(AP).lastnc14time = lastnc14time;
ResultsPlot(AP).t_raw = t(Frame_raw);
ResultsPlot(AP).MS2_raw = MS2_raw;
end

%% Save data into .mat structure
%Replace empty results with nans
for i = 1:N_bins
    fields = fieldnames(FitResults(i));
    for j = 1:length(fields)
        if isempty(FitResults(i).(fields{j}))
            FitResults(i).(fields{j}) = nan;
        end
    end
end

loc = [DropboxFolder,filesep,Prefix,'\'];
filename = 'MeanWholeCycleFitsNC14';
save([loc,filename,'.mat'],'FitResults','ResultsPlot','Prefix');

close(w);
disp(['Whole cycle fit analysis complete. Information stored in: ',loc]);

%Close the parallel pool if desired
if ~KeepPool
    poolobj = gcp('nocreate');
    delete(poolobj);
end
end

function [basal, t_on, t_dwell, t_active, tau, rate, residual, resnorm, CI] = FitWholeCycle(MS2,timepoints,t_interp)
%Fitting script to fit whole cycle to MS2 signal

%Anonymous fitting function

%Interpolate signal to data
fun = @(x)interp1(t_interp,MakeWholeCycle(x(1),x(2),x(3),x(4),x(5),x(6),t_interp),timepoints) - MS2;

%Fit the data
basal0 = 1;
t_on0 = timepoints(1);
t_dwell0 = 5;
t_active0 = 3;
tau0 = 10;
rate0 = 10;
x0 = [basal0, t_on0, t_dwell0, t_active0, tau0, rate0]; %Initial condition
%x_typ = [5, 3, 5, 7, 10, 10]; %typical values
lb = [0, 0, 0, -2, 0, 0]; %lower bound
ub = [100, 20, 20, 30, 30, 1e7]; %upper bound
options = optimoptions(@lsqnonlin,'Display','off','DiffMinChange',0.1);

[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(fun,x0,lb,ub,options); %Fit

%Extract fit parameters
basal = x(1);
t_on = x(2);
t_dwell = x(3);
t_active = x(4);
tau = x(5);
rate = x(6);
resnorm = resnorm ./ length(MS2); %Normalize residual norm by number of datapoints.

%Calculate confidence interval
CI=nlparci(x,residual,'jacobian',jacobian);

%Standard error in rate estimate


end

function signal = MakeWholeCycle(basal,t_on,t_dwell,t_active,tau,rate,t_interp)
%Function to generate a linear rise in signal to steady state
signal = zeros(size(t_interp)); %Output signal
dt = t_interp(2) - t_interp(1); %Time resolution
ind_shift = floor(t_dwell/dt); %Number of indices to shift for dwell time


%Generate loading rate
rate_cycle = nan(size(t_interp)); %Loading rate for the nuclear cycle
rate_cycle(t_interp < t_on) = 0;
rate_cycle(and(t_interp >= t_on, t_interp < (t_on + t_dwell + t_active))) = rate;
for i = 1:length(t_interp)
    if t_interp(i) >= (t_on + t_dwell + t_active)
        rate_cycle(i) = rate * exp(-(t_interp(i) - t_on - t_dwell - t_active)/tau);
    end
end

%Generate signal
for i = 2:length(t_interp)
    if t_interp(i) > t_on && t_interp(i) <= (t_on + t_dwell)
        signal(i) = signal(i-1) + dt * rate_cycle(i-1);
    elseif t_interp(i) > (t_on + t_dwell)
        if (i-1-ind_shift) < 1 %Check for zero indexing error
            ind_shift = ind_shift - 1;
        end
        signal(i) = signal(i-1) + (rate_cycle(i-1) - rate_cycle(i-1-ind_shift)) * dt;
    end
end

%Add basal fluorescence
signal(signal < basal) = basal;
end