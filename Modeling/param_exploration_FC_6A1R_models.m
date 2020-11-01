%% Parameter exploration
%% FC prediction at different regimes

% Predict the fold-change with sets of parameters ranging from 10^(-2) to
% 10^2 regime for Kb, Kr, wb, wbp, w_rp, w_br, w_brp, p (R_max).

% Flow : 
% 1) Use the Bcd and Run datasets for input TF gradients
% 2) Set the range of parameters,
% 3) Write down the functional form of FC for different models 
% 4) Calculate the fold-change with a specific set of parameters. Loop this
% over all combination of parameters.

% Then, it will give us three values, FC at 20%, 30%, and 40% of the embryo. 
% We will then plot these predictions in 2D plot with pairs of (20%, 40%), (20%, 30%), (30%, 40%) 


%% Define the filepaths
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo\FC_exploration';

%% Color module
% This is defining the line color
% We have 8 distinct datasets, with or without Runt protein.
% I think selecting 8 distinguishable color sets, then changing the
% brightness by either adding/subtracting white would be a better idea than
% selecting 16 different color sets.

colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.purple = [171,133,172]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

%colorDict.magenta = [208,109,171]/255;
%colorDict.lightBlue = [115,142,193]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.pink = [232,177,157]/255;
colorDict.thickpink = [132,27,69]/255;

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.darkgreen; colorDict.thickpink]; 

% For now, I'll add white (color+[1 1 1])/2 to make thinner color (for the
% Runt nulls)

%% Input data (Import the processed, time-averaged data)
TFData = load([FilePath, filesep, 'TFinput.mat']);
TFdata = TFData.TFinput;

Bcd = TFdata(:,1);
Run = TFdata(:,2);
RunNull = TFdata(:,3);

%% Import the FC data (from the compiledData.mat)
%% Define the parameter ranges

% 1) affinity (nM, AU)
K_B = logspace(-2,2,5);
K_R = logspace(-2,2,5);

% 2) interaction terms (exp(-eps/kBT)
W_b = logspace(0,2,5); % w_b
W_bp = logspace(0,2,5); % w_bp
% 3) interaction terms (repression)
W_br = logspace(-3,0,5); % w_br
W_rp = logspace(-3,0,5); % w_rp
W_brp = logspace(-3,0,5); % w_brp



% 3) RNAP/K_p
P = logspace(-3,-1,5); % p

% Initialize the fold-change matrix at 20%, 30%, and 40%
FC_20 = nan(length(K_B),length(K_R),length(W_b),length(W_bp),length(W_br),length(W_rp),length(W_brp),length(P));
FC_30 = nan(length(K_B),length(K_R),length(W_b),length(W_bp),length(W_br),length(W_rp),length(W_brp),length(P));
FC_40 = nan(length(K_B),length(K_R),length(W_b),length(W_bp),length(W_br),length(W_rp),length(W_brp),length(P));

APbin1 = 20/2.5+1;
APbin2 = 30/2.5+1;
APbin3 = 40/2.5+1;

% Loop through all varialbes to calculate the fold-change
for i= 1:length(K_B)
    for j= 1:length(K_R)
        for k = 1:length(W_b)
            for l= 1:length(W_bp)
                for m= 1:length(W_br)
                    for n= 1:length(W_rp)
                        for h= 1:length(W_brp)
                            for g= 1:length(P)
                                clear params
                                params = [K_B(i), K_R(j), W_b(k), W_bp(l), W_br(m), W_rp(n), W_brp(h),P(g)];
                                FC_20(i,j,k,l,m,n,h,g) = model_FC_6A1R_combination_all_V1(Bcd(APbin1), Run(APbin1), params);
                                FC_30(i,j,k,l,m,n,h,g) = model_FC_6A1R_combination_all_V1(Bcd(APbin2), Run(APbin2), params);
                                FC_40(i,j,k,l,m,n,h,g) = model_FC_6A1R_combination_all_V1(Bcd(APbin3), Run(APbin3), params);
                            end
                        end
                    end
                end
            end
        end
    end
end

%% plot the FC at different regions in 2D space
%% Plot FC at 20% versus 40% of the embryo

% step1. flatten the matrix 
FC_ant = reshape(FC_20,[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(W_rp)*length(W_brp)*length(P)]);
FC_mid = reshape(FC_30,[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(W_rp)*length(W_brp)*length(P)]);
FC_post = reshape(FC_40,[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(W_rp)*length(W_brp)*length(P)]);

%% Individual mode of repression : For each, we just set the interaction term, omega to be 1.
% First, competition : set w_rp = 1, w_brp = 1

FC_ant_compete = reshape(FC_20(:,:,:,:,:,end,end,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);
FC_mid_compete = reshape(FC_30(:,:,:,:,:,end,end,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);
FC_post_compete = reshape(FC_40(:,:,:,:,:,end,end,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);

% Second, direct repression : set w_br =1, w_brp = 1
FC_ant_direct = reshape(FC_20(:,:,:,:,end,:,end,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);
FC_mid_direct = reshape(FC_30(:,:,:,:,end,:,end,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);
FC_post_direct = reshape(FC_40(:,:,:,:,end,:,end,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);

% Third, quenching : set w_br = 1, w_rp = 1;
FC_ant_quench = reshape(FC_20(:,:,:,:,end,end,:,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);
FC_mid_quench = reshape(FC_30(:,:,:,:,end,end,:,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);
FC_post_quench = reshape(FC_40(:,:,:,:,end,end,:,:),[1,length(K_B)*length(K_R)*length(W_b)*length(W_bp)*length(W_br)*length(P)]);

%% generate the scatter plot of all possible FC predictions from each model, as well as real data points.
hold on
h(1) = scatter(FC_ant, FC_post, 'MarkerFaceColor', ColorChoice(1,:), 'MarkerEdgeColor', ColorChoice(1,:),'MarkerFaceAlpha',0.7)
h(2) = scatter(FC_ant_compete, FC_post_compete, 'MarkerFaceColor', ColorChoice(3,:), 'MarkerEdgeColor', ColorChoice(3,:),'MarkerFaceAlpha',0.7)
h(3) = scatter(FC_ant_direct, FC_post_direct, 'MarkerFaceColor', ColorChoice(4,:), 'MarkerEdgeColor', ColorChoice(4,:),'MarkerFaceAlpha',0.7)
h(4) = scatter(FC_ant_quench, FC_post_quench, 'MarkerFaceColor', ColorChoice(5,:), 'MarkerEdgeColor', ColorChoice(5,:),'MarkerFaceAlpha',0.7)

h(5) = errorbarxy(FC_001(APbin1), FC_001(APbin3),...
            FC_SEM_001(APbin1), FC_SEM_001(APbin3))
h(6) = errorbarxy(FC_010(APbin1), FC_010(APbin3),...
            FC_SEM_010(APbin1), FC_SEM_010(APbin3))
h(7) = errorbarxy(FC_100(APbin1), FC_100(16),...
            FC_SEM_100(APbin1), FC_SEM_100(16))

xlabel('FC-anterior')
ylabel('FC-posterior')

%% find the boundary of the scatter plots using the "boundary" function
% k = boundary(x,y) : k is the index of the boundary points, out of x-y
% vector. Note that x and y are column vectors, so we need transposition.

% The boundary function takes really long to calculate the boundary when
% the vectors have 10^8 elements. So, we'll sample the boundary with
% portion of the vectors, then sample again among those sampled-boundaries.
bound_temp = [];
for i=1:100
    % random sampling of 10^8 elements multiple times. We need to show that
    % the sampling number is enough to sample the whole space.
    index = randperm(10^8,10^4);
    bound_temp = unique([bound_temp ; boundary(FC_ant(index)',FC_post(index)')]);
end

% bound_temp = unique(bound_temp);

bound_combination = boundary(FC_ant(bound_temp)', FC_post(bound_temp)');

scatter(FC_ant(bound_temp)', FC_post(bound_temp)')

%%
bound_combination = boundary(FC_ant', FC_post');

bound_compete = boundary(FC_ant_compete', FC_post_compete');
bound_direct = boundary(FC_ant_direct', FC_post_direct');
bound_quench = boundary(FC_ant_quench', FC_post_quench');

%%
hold on
plot(FC_ant(bound_combination), FC_post(bound_combination),'color', ColorChoice(1,:),'LineWidth',1.5);
plot(FC_ant_compete(bound_compete), FC_post_compete(bound_compete), 'color', ColorChoice(2,:),'LineWidth',1.5)
plot(FC_ant_direct(bound_direct), FC_post_direct(bound_direct), 'color', ColorChoice(4,:),'LineWidth',1.5)
plot(FC_ant_quench(bound_quench), FC_post_quench(bound_quench), 'color', ColorChoice(6,:),'LineWidth',1.5)

errorbarxy(FC_001(APbin1), FC_001(APbin3),...
            FC_SEM_001(APbin1), FC_SEM_001(APbin3))
errorbarxy(FC_010(APbin1), FC_010(APbin3),...
            FC_SEM_010(APbin1), FC_SEM_010(APbin3))
errorbarxy(FC_100(APbin1), FC_100(16),...
            FC_SEM_100(APbin1), FC_SEM_100(16))

xlabel('FC-anterior')
ylabel('FC-posterior')

xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])

legend('combination','competition','diret','quenching','Location','NorthWest')

StandardFigure(gcf,gca)

% save the plot
saveas(gcf,[FigPath,filesep,'params_exploration_ant_post','.pdf']);

%% find the 3D-boundary

bound_combination = boundary(FC_ant', FC_mid', FC_post');

bound_compete = boundary(FC_ant_compete', FC_mid_compete', FC_post_compete');
bound_direct = boundary(FC_ant_direct', FC_mid_direct', FC_post_direct');
bound_quench = boundary(FC_ant_quench', FC_mid_quench', FC_post_quench');

%%
hold on
trisurf(bound_combination, FC_ant, FC_mid, FC_post,'Facecolor', ColorChoice(1,:),'FaceAlpha',0.5);
trisurf(bound_compete, FC_ant_compete, FC_mid_compete, FC_post_compete, 'Facecolor', ColorChoice(2,:),'FaceAlpha',0.5)
trisurf(bound_direct, FC_ant_direct, FC_mid_direct, FC_post_direct, 'Facecolor', ColorChoice(4,:),'FaceAlpha',0.5)
trisurf(bound_quench, FC_ant_quench, FC_mid_quench, FC_post_quench, 'Facecolor', ColorChoice(6,:),'FaceAlpha',0.5)

 view(-40,40)

plot3(FC_001(APbin1), FC_001(APbin2), FC_001(APbin3),'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
plot3(FC_010(APbin1), FC_010(APbin2), FC_010(APbin3),'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
plot3(FC_100(APbin1), FC_100(APbin2), FC_100(16),'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')

xlabel('FC-anterior')
ylabel('FC-mid')
zlabel('FC-posterior')

xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2])
ylim([0 1.2])
yticks([0 0.2 0.4 0.6 0.8 1 1.2])
zlim([0 1.2])
zticks([0 0.2 0.4 0.6 0.8 1 1.2])

legend('combination','competition','diret','quenching','Location','NorthWest')

box on
grid on

StandardFigure(gcf,gca)

% save the plot
saveas(gcf,[FigPath,filesep,'params_exploration_3D','.pdf']);