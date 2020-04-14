function [AccumulatedmRNA, AccumulatedmRNA_SD, AccumulatedmRNA_SE,...
            NC12,NC13,NC14,...
            AccumulatedmRNA_individual, AccumulatedmRNA_SD_individual] = accumulatemRNA(DataType, MinParticles, varargin)
% This function is for calculating the accumulated mRNA,
% using Processed datasets from AverageDatasets.m
% So, this script will
% 1) grabs MeanVectorAP, SDVectorAP, and NParticlesAP from individual
% embryos, then calculate the accumulated mRNA at each time point.
% 2) calculates the mean, SEM of accumulated mRNA over multiple embryos
%
% Input : 
% (1) DataType : Name of the dataset, after the AverageDatasets.m
% (2) 
% Option : 
% (1) mRNA half-life
% (2) ONnuclei : whether it's over ALL nuclei or ON nuclei
%  For ON Nuclei, 

%% initialize the variables (options)
ONnuclei = 0; % whether it's accumulated for only ON nuclei

%varargin = varargin{1};

% Checking Varargin 
if isempty(varargin)
    % default
elseif strcmpi(varargin{1},'ONnuclei')
    % Accumulate over ON nuclei only
    ONnuclei = 1;
end
%% Load the processed data
% Define the file path. This should be fed into the function as an optional
% input. Edit later.
% FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData'; 
% Dropbox filePath
DropboxPath = 'S:/YangJoon/Dropbox';
FilePath = [DropboxPath,filesep,'OpposingGradient/OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];
Data = load([FilePath, filesep, DataType,'.mat']);

%% Extract useful fields
% number of embryos used
numEmbryos = length(Data.MeanVectorAP_individual(1,1,:));

ElapsedTime = Data.ElapsedTime;
NC12 = Data.nc12;
NC13 = Data.nc13;
NC14 = Data.nc14;

MeanVectorAP = Data.MeanVectorAP_individual;
SDVectorAP = Data.SDVectorAP_individual;
NParticlesAP = Data.NParticlesAP_individual;

%% Convert NaNs to Zeros
MeanVectorAP(isnan(MeanVectorAP)) = 0;
SDVectorAP(isnan(SDVectorAP)) = 0;
NParticlesAP(isnan(NParticlesAP)) = 0;

%% Calculate the accumulated mRNA

NFrames = length(ElapsedTime);
nAPbins = 41; % This should be done better.

AccumulatedmRNA_individual = zeros(NFrames,nAPbins,numEmbryos);
AccumulatedmRNA_SD_individual =  zeros(NFrames,nAPbins,numEmbryos);

% Loop over individual embryos
for k=1:numEmbryos
    % Loop over all AP bins 
    % Here, I need to consider the APbinArea, 
    % like small AP bins should be ignored, since it'll only add bias.
    for i=1:nAPbins
        for j=2:length(ElapsedTime)
            % filter for the minimum number of particles in each AP bin.
            if NParticlesAP(1:j,i,k)<MinParticles
                NParticlesAP(1:j,i,k) = 0;
            end
                
            % OnNuclei check
            if ONnuclei==0
                AccumulatedmRNA_individual(j,i,k) =...
                    trapz(ElapsedTime(1:j),MeanVectorAP(1:j,i,k).*NParticlesAP(1:j,i,k)); 
                AccumulatedmRNA_SD_individual(j,i,k) =...
                    sqrt(trapz(ElapsedTime(1:j),SDVectorAP(1:j,i,k).^2 .*NParticlesAP(1:j,i,k)));
            elseif ONnuclei ==1
                AccumulatedmRNA_individual(j,i,k) =...
                    trapz(ElapsedTime(1:j),MeanVectorAP(1:j,i,k)); 
                AccumulatedmRNA_SD_individual(j,i,k) =...
                    sqrt(trapz(ElapsedTime(1:j),SDVectorAP(1:j,i,k).^2 ));
                
        end
    end
end

%% Convert Zeros to NaNs back
AccumulatedmRNA_individual(AccumulatedmRNA_individual==0) = nan;
AccumulatedmRNA_SD_individual(AccumulatedmRNA_SD_individual==0) = nan;

%% Plot to check
% tPoint = 1;
% 
% for tPoint = 1:length(ElapsedTime)
%     clf
%     hold on
%     for i=1:numEmbryos
%         errorbar(0:0.025:1, AccumulatedmRNA_individual(tPoint, :,i), AccumulatedmRNA_SD_individual(tPoint, :,i))
%     end
%     pause
% end
%% Average over multiple embryos
AccumulatedmRNA = nanmean(AccumulatedmRNA_individual,3);
AccumulatedmRNA_SD = nanstd(AccumulatedmRNA_individual,0,3);
AccumulatedmRNA_SE = AccumulatedmRNA_SD./sqrt(numEmbryos);

%% Plot to check, individual vs mean
% for tPoint = 1:length(ElapsedTime)
%     clf
%     hold on
%     for i=1:numEmbryos
%         errorbar(0:0.025:1, AccumulatedmRNA_individual(tPoint, :,i), AccumulatedmRNA_SD_individual(tPoint, :,i))
%     end
%     errorbar(0:0.025:1, AccumulatedmRNA(tPoint,:), AccumulatedmRNA_SE(tPoint,:))
%     pause
% end

%% Save fields
% DropboxPath = 'S:/YangJoon/Dropbox';
% FilePath = [DropboxPath,filesep,'OpposingGradient/OpposingGradients_ProcessedData'];
% savedVariables = {};
% savedVariables = [savedVariables, 'AccumulatedmRNA', 'AccumulatedmRNA_SD',...
%                 'AccumulatedmRNA_SE','NC12','NC13','NC14',...
%                 'AccumulatedmRNA_individual', 'AccumulatedmRNA_SD_individual'];
% 
% save([FilePath,filesep,'AccumulatedmRNA_',DataType,'.mat'],...
%     savedVariables{:},'-v7.3');

end