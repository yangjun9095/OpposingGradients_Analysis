function main20_Check_MCP_Dosage_Datasets(DataType)
% This script is for checking the MCP dosage across multiple constructs.
% The logic is below : 
% For datasets with MS2 spots, we can use the offset values used for spot
% intensity calculation as a proxy for dosage of MCP-FP.
% CompiledParticles.MeanOffsetVector
% Note. MeanOffsetVector has dimension of "1 x (frames)"
% Last updated : Nov.2019
%% Define the datasets that we want to compare
% Make a string with DataTypes
DataString = {'r0','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3'};
%% Load the datasets
% Data = LoadMS2Sets(DataType);
%% Time average the MeanOffsetVector to compare across embryos within a dataset
% % Caveats : There's some fluctuation in MCP-FP nuclear level over time, but
% % we'll just average over cycles to get some rough number we can compare
% % across datasets.
% 
% for i=1:length(Data)
%     Offset_mean(i) = nanmean(Data(i).MeanOffsetVector);
%     Offset_STD(i) = nanstd(Data(i).MeanOffsetVector);
%     Offset_SEM(i) = Offset_STD(i)./sqrt(length(Data(i).MeanOffsetVector));
% end
%     
% errorbar(1:4, Offset_mean, Offset_STD)

%% Averaged Offset values across DataTypes
% Go through each DataType 
for i=1:length(DataString)
    Data = LoadMS2Sets(DataString{i});
    % go through each embryo within the same DataType
    for j=1:length(Data)
        Offset_mean_ind(j) = nanmean(Data(j).MeanOffsetVector);
    end
    Offset_mean(i) = nanmean(Offset_mean_ind);
    Offset_STD(i) = nanstd(Offset_mean_ind);
    Offset_SEM(i) = Offset_STD(i)./sqrt(length(Data));
end

errorbar(1:length(DataString), Offset_mean, Offset_SEM)
ylim([0 12])
set(gca,'XtickLabel',DataString)
title('MCP-GFP dosage check')
xlabel('Constructs')
ylabel('MCP-GFP dosage (Offset) (AU)')
StandardFigure(gcf,gca)

% Save plots
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\MCP-GFP_DosageCheck';
% 
% saveas(gcf,[FigPath,filesep,'MCP-GFP_Dosage_Check_AllConstructs','.tif'])
% saveas(gcf,[FigPath,filesep,'MCP-GFP_Dosage_Check_AllConstructs','.pdf'])

%% Save the processed MCP-GFP dosage data
savedVariables = {};
savedVariables = [savedVariables,'DataString',...
                                    'Offset_mean',...
                                    'Offset_SEM'];
                                
% Save the variables into .mat file.
FilePath = FigPath; 
save([FilePath,filesep,'MCP-GFP_dosage_AllConstructs.mat'],...
    savedVariables{:},'-v7.3');

end