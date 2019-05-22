function main_11_Accumulated_mRNA_sex_dependence(DataType, varargin)
%% DESCRIPTION
% This script is to compare the integrated mRNA profiles (over AP, at
% different NC ranges, etc.) between males and females of the same
% synthetic enhancers.

% The plan is to use this script to compare 
% 1) r0, show that they show equivalent expression (quantitatively)
% 2) r3 (r1, r2, and r3' as well) to show that how much of sex-dependence
% is present in here.

%% Load datasets
% This assumes that the DataType is the name of the constructs in the
% DataStatus.xlsx tab, for example, r3, and also has either '-male', or
% '-female' after that DataType.

males = LoadMS2Sets([DataType,'-male']);
females = LoadMS2Sets([DataType,'-female']);

%% Calculate the total mRNA (Accumulated mRNA, or IntegratedmRNA)
% There are multiple ways to calculate this right now.
% 1) HG's IntegratemRNA.m 
% 2) YJK's AverageDatasets.m
% 3) 

% For now, let's use HG's IntegratemRNA script
[TotalProd_male,TotalProdError_male,TotalProdN_male,...
    MeanTotalProd_male,SDTotalProd_male,SETotalProd_male] = ...
                        IntegratemRNA(males,2,1)
                    
[TotalProd_female,TotalProdError_female,TotalProdN_female,...
    MeanTotalProd_female,SDTotalProd_female,SETotalProd_female] = ...
                        IntegratemRNA(females,2,1)
                    
hold on                    
for i=1:length(males)
    errorbar(0:0.025:1, TotalProd_male(i,:,13), TotalProdError_male(i,:,13),'b')
    pause
end

for i=1:length(females)
    errorbar(0:0.025:1, TotalProd_female(i,:,13), TotalProdError_female(i,:,13),'r')
    pause
end

ylim([0 8000])
% errorbar(0:0.025:1,MeanTotalProd_male(:,13),SETotalProd_male(:,13))
% errorbar(0:0.025:1,MeanTotalProd_female(:,13),SETotalProd_female(:,13))
end