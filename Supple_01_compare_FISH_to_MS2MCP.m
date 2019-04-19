function Supple_01_compare_FISH_to_MS2MCP
% DESCRIPTION
% This script is for comparison of FISH and MS2-MCP in terms of accumulated
% cytoplasmic mRNA pattern. 
% As the first step, we will compare Jeehae's FISH data for P2P, with our
% P2P-MS2-MCP data. Should we try MS2.V5? I think it's good to use MS2.V5.

%% Define the directory
%gives the location of these 5 quantities.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
%% Load the P2P-lacZ FISH datasets
FISHPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\hbP2P-LacZ-FISH_LineTraces_Jeehae\wtHbP2-LacZ in situ traces for Yangjoon Kim';
D = dir(FISHPath);
for i=3:length(D)
    FISH_data(i-2,:) = load(['E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\hbP2P-LacZ-FISH_LineTraces_Jeehae\wtHbP2-LacZ in situ traces for Yangjoon Kim\',D(i).name]);
end

%% Plot the data over AP
hold on
for i=1:length(FISH_data(:,1))
    plot(0.005:0.01:0.995, FISH_data(i,:))
    %pause
end

title ('FISH pattern over AP')
xlabel('AP axis (EL)')
ylabel('Normalized intensity (AU)')

%% Average the FISH data over embryos
[numEmbryos_FISH, numAPbins_FISH] = size(FISH_data);
AP_FISH = 0.005:0.01:0.995;

Mean_FISHData = nanmean(FISH_data);
SEM_FISHData = nanstd(FISH_data)./sqrt(numEmbryos_FISH);

errorbar(AP_FISH, Mean_FISHData, SEM_FISHData)
%% Load the P2P-MS2.V5-lacZ datasets
P2P_MS2V5_data = LoadMS2Sets('P2P-MS2v5-lacZ-36uW')
P2P_MS2V1_data = LoadMS2Sets('P2P-MS2-lacZ-36uW')
%% Calculate the total mRNA (accumulated mRNA) for P2P-MS2V5 datasets using IntegratemRNA.m
% [TotalProd,TotalProdError,TotalProdN,...
%     MeanTotalProd,SDTotalProd,SETotalProd]=IntegratemRNA(Data,MinParticles,MinEmbryos,varargin)
[TotalProd,TotalProdError,TotalProdN,...
    MeanTotalProd,SDTotalProd,SETotalProd]=IntegratemRNA(P2P_MS2V5_data,2,2)

%% plot the Mean total Prod (this is total accumulated mRNA over specific NC, divided by the total number of nuclei (ON + Off)
% Since these MeanTotalProd is divided with total # of nuclei, I will
% multiply with the scaling factor to account for the different numbers of
% nuclei in different cycles.

% Scale factor : Assuming that all nuclei undergo doubling at each cycle.
%scale_NC11 = 0.125;
scale_NC12 = 0.25;
scale_NC13 = 0.5;
scale_NC14 = 1;

hold on
%errorbar(0:0.025:1, MeanTotalProd(:,11)*scale_NC11, SETotalProd(:,11)*scale_NC11)
errorbar(0:0.025:1, MeanTotalProd(:,12)*scale_NC12, SETotalProd(:,12)*scale_NC12)
errorbar(0:0.025:1, MeanTotalProd(:,13)*scale_NC13, SETotalProd(:,13)*scale_NC13)
errorbar(0:0.025:1, MeanTotalProd(:,14)*scale_NC14, SETotalProd(:,14)*scale_NC14 )

title('Accumulated mRNA over AP')
xlabel('AP (EL)')
ylabel('Accumulated mRNA (AU)')
legend('NC12','NC13','NC14')

StandardFigure(gcf,gca)
%% Try to calculate the total mRNA over NC13 and NC14

% Scale factor : Assuming that all nuclei undergo doubling at each cycle.
%scale_NC11 = 0.125;
scale_NC12 = 0.25;
scale_NC13 = 0.5;
scale_NC14 = 0;

MeanTotalProd_NC12NC14 = MeanTotalProd(:,12)*scale_NC12 +...
                            MeanTotalProd(:,13)*scale_NC13 + MeanTotalProd(:,14) * scale_NC14; %MeanTotalProd(:,11)*scale_NC11 + 
SETotalProd_NC12NC14 = sqrt (SETotalProd(:,12).^2*scale_NC12 + ...
                                SETotalProd(:,13).^2*scale_NC13 + SETotalProd(:,14)* scale_NC14);%SETotalProd(:,11).^2*scale_NC11 + 

% Subtract the background (I should do this better)
MeanTotalProd_NC12NC14 = MeanTotalProd_NC12NC14 - min(MeanTotalProd_NC12NC14);
                            
MeanTotalProd_NC12NC14_Normalized = MeanTotalProd_NC12NC14 ./ max(MeanTotalProd_NC12NC14);
SETotalProd_NC12NC14_Normalized = SETotalProd_NC12NC14./max(MeanTotalProd_NC12NC14);
%% Normalize the Mean Total Prod for NC13 and NC14
MeanTotalProd_NC13_Normalized = MeanTotalProd(:,13) ./ max(MeanTotalProd(:,13));
SETotalProd_NC13_Normalized = SETotalProd(:,13)./ max(MeanTotalProd(:,13));

MeanTotalProd_NC14_Normalized = MeanTotalProd(:,14) ./ max(MeanTotalProd(:,14));
SETotalProd_NC14_Normalized = SETotalProd(:,14) ./ max(MeanTotalProd(:,14));

%% plot both MS2 and FISH data together
ComparisonFigure = figure;
hold on
errorbar(AP_FISH, Mean_FISHData, SEM_FISHData)
errorbar(0:0.025:1, MeanTotalProd_NC12NC14_Normalized*1.4,SETotalProd_NC12NC14_Normalized)

title('Comparison of MS2-MCP and FISH - P2P')
xlabel('AP axis (EL)')
ylabel('Accumulated mRNA (normalized)')
legend('FISH-Jeehae','MS2V5')
StandardFigure(gcf,gca)

FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\DynamicsResults\hbP2P-LacZ-FISH_LineTraces_Jeehae\wtHbP2-LacZ in situ traces for Yangjoon Kim';
%  saveas(ComparisonFigure,[FigPath,filesep,'Comparison_AccumulatedmRNA_FISH_MS2_from_NC12_toNC14_BackgroundSubtracted' , '.tif']); 
%  saveas(ComparisonFigure,[FigPath,filesep,'Comparison_AccumulatedmRNA_FISH_MS2_from_NC12_toNC14_BackgroundSubtracted' , '.pdf']); 

%% Save the variables
save(['E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData',filesep,'Comparison_FISH_MS2V5-P2P-lacZ.mat'],...
                            'FISH_data','AP_FISH','Mean_FISHData','SEM_FISHData',...
                            'P2P_MS2V5_data','MeanTotalProd_NC12NC14','MeanTotalProd_NC12NC14_Normalized',...
                            'SETotalProd_NC12NC14', 'SETotalProd_NC12NC14_Normalized','-v7.3')
end