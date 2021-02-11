%% process data for the MCMC function
function process_inputTF_outputRates_1RunSite_forMCMC
%% Description

%% Load the datasets
%% Import data for the MCMC inference
% From the "preprocess_data_for_MCMC.m" script
% xdata(TFinputs) and ydata(initial rate), note that we will do
% a simultaneous fitting for the Runt WT and Runt nulls.

% We need another separate script to process the data for inputs in this
% script. : This is now done in the "preprocess_data_for_MCMC.m" script.

preprocessedData = load([FilePath, filesep, 'PreProcessedData_ForMCMC.mat']);
data = preprocessedData.data;

% Load the input TF data
load([FilePath, filesep, 'TFinput.mat'])

Bicoid = TFinput(:,1);
Runt = TFinput(:,2);
RuntNull = TFinput(:,3);

%% Pull the construct of our interest (for the parameter inference)
% initialize the MCMCdata structure 
% MCMCdata = cell;

% initialize the counter
k=1; 

for construct = [2,5,6]
%     % Choose a construct 
%     construct = 2;
    % Pick the dataset from the data.mat
    Data = data(construct);

    % MCMC analysis on the initial slope (averaged over embryos)
    % initialize the AP bins
    APaxis = Data.APbins;

    %Truncate APbins to user-specified range (input was optional argument in
    %this function.
    NoNaN_index_null = ~isnan(Data.Rate_null);
    NoNaN_index_WT = ~isnan(Data.Rate_WT);
    % calculate the AP bins that are not NaNs in both WT and Null datasets
    NoNaN_index = NoNaN_index_null.*NoNaN_index_WT;

    NoNaNindices = find(NoNaN_index);

    % Range that is set as an initial guess. We will get a common set of
    % APbins that does not have NaN values in these AP bins.
    APbin_start = 20/2.5 + 1;
    APbin_end = 45/2.5 + 1;

    APbinRange = (APbin_start:APbin_end)';

    % find the common elements of AP bins between Not-NaNs, and also the
    % pre-set range (20-45%)
    APbins_fit = intersect(NoNaNindices, APbinRange);

    % initialize the initial rate (slope)
    Rate_WT = Data.Rate_WT;
    Rate_null = Data.Rate_null;

    % Truncate the vectors using the range of AP bins
    APbins = APaxis(APbins_fit);
    Rate_WT = Rate_WT(APbins_fit);
    Rate_null = Rate_null(APbins_fit);

    Bcd = Bicoid(APbins_fit);
    Run = Runt(APbins_fit);
    RunNull = RuntNull(APbins_fit);

    % Decide whether we want to fit the Runt null/WT data together or not.
    % depending on this, we will set the xdata and ydata for the fitting.
%     MCMCdata = struct;
    MCMCdata{k}.APdata = [APbins' ; APbins'];
    MCMCdata{k}.ydata = [Rate_null ; Rate_WT];
    
    % define the TFtemp to parse the TF input info (both Bcd and Run)
    % input TF
%     MCMCdata{k}.Bcd = [Bcd ; Bcd];
%     MCMCdata{k}.Run = [RunNull ; Run];
    clear TFtemp
    TFtemp(:,1) = [Bcd; Bcd];
    TFtemp(:,2) = [RunNull; Run];
    MCMCdata{k}.xdata =  TFtemp;

    MCMCdata{k}.R_max = max([Rate_null ; Rate_WT]);
    MCMCdata{k}.R_min = min([Rate_null ; Rate_WT]);
    
    k=k+1;
end

%% save the MCMCdata structure for repeated use of different types of models
FilePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\MCMC_HillV3';
save([FilePath, filesep,'MCMCdata_1RunSite.mat'],...
        'MCMCdata')
end