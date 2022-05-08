function PATSER_02_generate_Runt_sites_variants
%% Description
% generate plots of PATSER scores for a given enhancer for given sets of
% TFs. We will utilize the mapPATSERResults.m script for this.

%% Directories
PWMpath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PATSER\pwm_Park&DePace_eLife_2020';
DataPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PATSER';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PATSER\PATSER_figures';

%% Runt site variants and their PATSER scores
% Refer to the Google Docs on which sequences were run for PATSER, and the
% justifications.

index_Constructs = [1 2 3 4 5 6 7 8 9 10];
name_constructs = {'PWM','','1','2','3','4','5','6','7','8'}
PATSER_scores = [6.51 4.62 3.02 4.15 3.60 5.59 1.39 3.42 3.55 3.04];

bar(index_Constructs, PATSER_scores,0.5)
xticklabels(name_constructs)
xlim([0 11])
ylim([2 7])
yticks([2 3 4 5 6 7])
% PATSER score from the wt hb P2 (which was NOT modulated by Runt protein)
% we choose the higest one as the lower limit for the real Runt binding
% site.
ID_var = 'hbP2_Tao';
PWM_code = 'Lab_repo';
Run_pos = [19 55 67 93 170 177 184 216];
Run_width = 9;
Run_PATSER_score = [3.20 3.38 3.27 3.58 3.88 3.78 3.93 3.22];
Run_PATSER_score_falsepositive = max(Run_PATSER_score);

yline(Run_PATSER_score_falsepositive,'--')

box on
StandardFigure(gcf,gca)

% save the plot
saveas(gcf,[FigPath, filesep, 'PATSER_','Runt_sites_variants_V1','.pdf'])

end