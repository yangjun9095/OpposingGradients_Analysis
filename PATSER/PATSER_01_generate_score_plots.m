%% generate plots of PATSER score results for a given enhancer for each individual TF
function PATSER_01_generate_score_plots(enhancer, TFs)
%% Description
% generate plots of PATSER scores for a given enhancer for given sets of
% TFs. We will utilize the mapPATSERResults.m script for this.

%% Temporary description
% For now, we will hard-code the PATSER result into list/vector of position, width, and scores for individual TFs from individual PWMs.

TF_pos = [];
TF_width = [];
TF_PATSER_score = [];

%% 1) FlyReg
Bcd_pos = [41 83 99 152 182 257];
Bcd_width = 8;
Bcd_PATSER_score = [3.65 3.94 3.74 3.88 3.46 4.37];

%% 2) Cell,2008
Bcd_pos = [41 83 99 152 182 257];
Bcd_width = 8;
Bcd_PATSER_score = [3.65 3.94 3.74 3.88 3.46 4.37];

%% 3) new5 NAR
Bcd_pos = [12 40 109 117 151 164 174 181 195 256];
Bcd_width = 8;
Bcd_PATSER_score = [3.62 4.62 3.61 3.56 4.53 3.85 3.42 3.78 3.49 4.53];

%% 4) SOLEXA


%% generate plots
enhancer_length = 294;
% Bcd_width_perc = Bcd_width/enhancer_length;

bar(Bcd_pos, Bcd_PATSER_score, 0.5)
 
xlim([0 294])
xticks([0 50 100 150 200 250 294])
ylim([3 5])
yticks([3 3.5 4 4.5 5])
box on

StandardFigure(gcf,gca)

% save the plot

end