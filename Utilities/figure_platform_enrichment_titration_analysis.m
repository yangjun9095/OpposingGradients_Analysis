clear
close all
addpath('../utilities')

% set paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
project1x = 'Dl-Ven_snaBAC-mCh';
project2x = '2xDl-Ven_snaBAC-mCh';
[~, DataPath1x, FigureRoot] =   header_function(DropboxFolder, project1x);
[~, DataPath2x, ~] =   header_function(DropboxFolder, project2x);
FigPath = [FigureRoot 'dorsal_sna_enrichment_titration/'];
mkdir(FigPath)
% load data
load([DataPath1x 'nucleus_struct_protein.mat'])
nc_data_1x = nucleus_struct_protein;
load([DataPath2x 'nucleus_struct_protein.mat'])
nc_data_2x = nucleus_struct_protein;
clear nucleus_struct_protein;

% basic plot and data qc params 
DistLim = 0.8; % min distance from edge permitted (um)
nBoots = 100;%00; % number of bootstrap samples to use for estimating SE

PixelSize = nc_data_1x(1).PixelSize;
zStep = nc_data_1x(1).zStep;
VoxelSize = PixelSize^2 * zStep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull useful vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from nucleus edge
dist1x_vec = [nc_data_1x.spot_edge_dist_vec];
dist1x_ft_vec = dist1x_vec >= DistLim;
dist2x_vec = [nc_data_2x.spot_edge_dist_vec];
dist2x_ft_vec = dist2x_vec >= DistLim;

% Dorsal at spot
spot_protein1x_vec = ([nc_data_1x.spot_protein_vec])*VoxelSize;
spot_protein1x_vec = spot_protein1x_vec(dist1x_ft_vec);
spot_protein2x_vec = [nc_data_2x.spot_protein_vec];
spot_protein2x_vec = spot_protein2x_vec(dist2x_ft_vec);

% absolute spot enrichment
delta_protein1x_vec = ([nc_data_1x.spot_protein_vec] - [nc_data_1x.edge_null_protein_vec])*VoxelSize;
delta_protein1x_vec = delta_protein1x_vec(dist1x_ft_vec);
delta_protein2x_vec = [nc_data_2x.spot_protein_vec] - [nc_data_2x.edge_null_protein_vec];
delta_protein2x_vec = delta_protein2x_vec(dist2x_ft_vec);

% average nuclear concentrations
mf_protein1x_vec = [nc_data_1x.mf_null_protein_vec]*VoxelSize;
mf_protein1x_vec = mf_protein1x_vec(dist1x_ft_vec);
mf_protein2x_vec = [nc_data_2x.mf_null_protein_vec];
mf_protein2x_vec = mf_protein2x_vec(dist2x_ft_vec);

% spot fluorescence
fluo1x_vec = [nc_data_1x.fluo];
fluo1x_vec = fluo1x_vec(dist1x_ft_vec);
fluo2x_vec = [nc_data_2x.fluo];
fluo2x_vec = fluo2x_vec(dist2x_ft_vec);

% time
time1x_vec = [nc_data_1x.time];
time1x_vec = time1x_vec(dist1x_ft_vec);
time2x_vec = [nc_data_2x.time];
time2x_vec = time2x_vec(dist2x_ft_vec);

%% Simple histograms to compare nucleus concentrations for 1 and 2x
close all
n_bins = 100;
% get bounds
mf_lb = min([mf_protein1x_vec mf_protein2x_vec]);
mf_ub = max([mf_protein1x_vec mf_protein2x_vec]);
dt_lb = prctile([delta_protein1x_vec mf_protein2x_vec],1);
dt_ub = prctile([delta_protein1x_vec mf_protein2x_vec],99);
% generate bins
delta_bins = linspace(dt_lb,dt_ub,n_bins);
mf_bins = linspace(mf_lb,mf_ub,n_bins);

% make figure2
mf_hist = figure;
cmap1 = brewermap([],'Set2');
cmap2 = brewermap([],'Paired');
hold on
histogram(mf_protein1x_vec,mf_bins,'Normalization','probability','FaceColor','blue','EdgeAlpha',.2)
histogram(mf_protein2x_vec,mf_bins,'Normalization','probability','FaceColor','red','EdgeAlpha',.2)
xlabel('nuclear Dl concentration (au)')
ylabel('share')
legend('1x','2x')
grid on 
box on
set(gca,'Fontsize',14)
saveas(mf_hist,[FigPath 'nuclear_dorsal_1x_vs_2x.png'])
saveas(mf_hist,[FigPath 'nuclear_dorsal_1x_vs_2x.pdf'])


%% Look at enrichment as a function of nucleus concentration

% calculate average enrichment as a function of nuclear concentration
% combine sets for now
% concatenate 1x and 2x sets
mf_vec_full = [mf_protein1x_vec mf_protein2x_vec];
time_vec_full = [time1x_vec time2x_vec];
delta_vec_full = [delta_protein1x_vec delta_protein2x_vec];
fluo_vec_full = [fluo1x_vec fluo2x_vec];
id_vec_full = [repelem(1,numel(mf_protein1x_vec)) repelem(2,numel(mf_protein2x_vec))];

% generate protein bins
n_points = 25; % number of protein bins
mf_lb = prctile([mf_protein1x_vec mf_protein2x_vec],.5);
mf_ub = prctile([mf_protein1x_vec mf_protein2x_vec],99.5);
mf_plot_index = linspace(mf_lb,mf_ub,n_points);

% generate vectors to use for filtering observations
type_cell = {1,2,1:2};
type_names = {'1x only', '2x only', 'combined'};
n_types = numel(type_names);
% initialize boot arrays
mf_sigma = median(diff(mf_plot_index))/2;
delta_boot_array = NaN(nBoots,n_points,n_types);
fluo_boot_array = NaN(nBoots,n_points,n_types);
time_boot_array = NaN(nBoots,n_points,n_types);

for t = 1:n_types
    type_ft = ismember(id_vec_full,type_cell{t});    
    for n = 1:n_points
        mf_pt = mf_plot_index(n);
        delta_mf_vec = mf_pt-mf_vec_full;
        mf_weights = exp(-.5*((delta_mf_vec)/mf_sigma).^2); % gaussian weights   
        mf_weights(abs(delta_mf_vec)>2*mf_sigma) = 0;
        type_options = find(mf_weights>0&type_ft);
        if numel(type_options) > 50
            for b = 1:nBoots
                boot_indices = randsample(type_options,numel(type_options),true);        
                boot_weights = mf_weights(boot_indices);                 
                delta_boot_array(b,n,t) = nansum(delta_vec_full(boot_indices).*boot_weights) / nansum(boot_weights);      
                fluo_boot_array(b,n,t) = nansum(fluo_vec_full(boot_indices).*boot_weights) / nansum(boot_weights);                
                time_boot_array(b,n,t) = nansum(time_vec_full(boot_indices).*boot_weights) / nansum(boot_weights);                
            end
        end
    end
end

%% plot
% calculate mean and se
mean_time_mean = squeeze(nanmean(time_boot_array));
mean_time_se = squeeze(nanstd(time_boot_array));

mean_enrichment_mean = squeeze(nanmean(delta_boot_array));
mean_enrichment_se = squeeze(nanstd(delta_boot_array));

mean_fluo_mean = squeeze(nanmean(fluo_boot_array));
mean_fluo_se = squeeze(nanstd(fluo_boot_array));

%%

mean_titration = figure;
hold on
s = [];
for n = 1:n_types
    e = errorbar(mf_plot_index,mean_enrichment_mean(:,n),mean_enrichment_se(:,n),'--','Color','black','LineWidth',1);
    e.CapSize = 0;
    s = [s scatter(mf_plot_index,mean_enrichment_mean(:,n),'MarkerFaceColor',cmap1(n,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','black')];
end
ylabel('Dl enrichment at locus (au)')
xlabel('nuclear Dl concentration (au)')
legend(s,type_names{:},'Location','northwest');
grid on 
box on
set(gca,'Fontsize',14)
saveas(mean_titration,[FigPath 'avg_titration_trend.png'])
saveas(mean_titration,[FigPath 'avg_titration_trend.pdf'])


%% Fit polynomial
close all
nan_filter = ~isnan(delta_vec_full) & ~isnan(mf_vec_full);
p3 = polyfit(mf_vec_full(nan_filter), delta_vec_full(nan_filter), 3);
p2 = polyfit(mf_vec_full(nan_filter), delta_vec_full(nan_filter), 3);

p2_trend = polyval(p2,mf_plot_index);
p3_trend = polyval(p3,mf_plot_index);

fit_titration = figure;
hold on
pl2 = plot(mf_plot_index,p2_trend,'-','LineWidth',1.5);
e = errorbar(mf_plot_index,mean_enrichment_mean,mean_enrichment_se,'o','Color','black' ,'LineWidth',1);
e.CapSize = 0;
s = scatter(mf_plot_index,mean_enrichment_mean,'MarkerFaceColor',cmap1(2,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','black');

% pl3 = plot(mf_plot_index,p3_trend,'-d','LineWidth',1.5);
legend([s pl2], 'mean','2nd order polynomial','Location','northwest')
ylabel('Dl enrichment at locus (au)')
xlabel('nuclear Dl concentration (au)')
grid on 
box on
set(gca,'Fontsize',14)
saveas(fit_titration,[FigPath 'fit_titration_trend.png'])
saveas(fit_titration,[FigPath 'fit_titration_trend.pdf'])

%% Look at fluoresncence vs. dorsal
close all
% just fluorescence
fluo_titration = figure;
fluo_color = [161 133 161]/256;
hold on
e = errorbar(mf_plot_index,mean_fluo_mean,mean_fluo_se,'-','Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(mf_plot_index,mean_fluo_mean,'MarkerFaceColor',fluo_color,'MarkerEdgeAlpha',1,'MarkerEdgeColor','black');

% pl3 = plot(mf_plot_index,p3_trend,'-d','LineWidth',1.5);
% legend([s pl2], 'mean','2nd order polynomial','Location','northwest')
ylabel('{\it snail} fluorescence (au)')
xlabel('nuclear Dl concentration (au)')
grid on 
box on
set(gca,'Fontsize',14)
saveas(fluo_titration,[FigPath 'fluo_titration_trend.png'])
saveas(fluo_titration,[FigPath 'fluo_titration_trend.pdf'])


% combined
cb_titration = figure;
hold on

yyaxis left
e = errorbar(mf_plot_index,mean_enrichment_mean,mean_enrichment_se,'--','Color','black' ,'LineWidth',1);
e.CapSize = 0;
scatter(mf_plot_index,mean_enrichment_mean,'MarkerFaceColor',cmap1(2,:),'MarkerEdgeAlpha',1,'MarkerEdgeColor','black');
ylabel('Dl enrichment at locus (au)')
ax = gca;
ax.YColor = 'black';%cmap1(2,:);

yyaxis right
e = errorbar(mf_plot_index,mean_fluo_mean,mean_fluo_se,'-','Color','black','LineWidth',1);
e.CapSize = 0;
scatter(mf_plot_index,mean_fluo_mean,'MarkerFaceColor',fluo_color,'MarkerEdgeAlpha',1,'MarkerEdgeColor','black');
ylabel('{\it snail} fluorescence (au)')
ax = gca;
ax.YColor = 'black';%'black';
xlabel('nuclear Dl concentration (au)')
grid on 
box on
set(gca,'Fontsize',14)
saveas(cb_titration,[FigPath 'combined_titration_trend.png'])
saveas(cb_titration,[FigPath 'combined_titration_trend.pdf'])

%% now plot fluo vs enrichment
close all
fluo_vs_delta = figure;
hold on
e1 = errorbar(mean_enrichment_mean,mean_fluo_mean,mean_fluo_se,'.','Color',fluo_color,'LineWidth',1.5);
e1.CapSize = 0;
e2 = errorbar(mean_enrichment_mean,mean_fluo_mean,-mean_enrichment_se*.001,...
    mean_enrichment_se*.001,-mean_enrichment_se,mean_enrichment_se,'.','Color',cmap1(2,:),'LineWidth',1.5);
e2.CapSize = 0;
scatter(mean_enrichment_mean,mean_fluo_mean,10,'MarkerFaceColor','black','MarkerEdgeAlpha',0,'MarkerFaceAlpha',1)
ylabel('{\it snail} fluorescence (au)')
xlabel('Dl enrichment at locus (au)')
grid on 
box on
set(gca,'Fontsize',14)
saveas(fluo_vs_delta,[FigPath 'fluo_vs_delta.png'])
saveas(fluo_vs_delta,[FigPath 'fluo_vs_delta.pdf'])