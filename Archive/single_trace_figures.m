% Script to generate figures comparing MS2 fitting methodologies
clear 
close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
% DropboxFolder = 'S:\Nick\Dropbox\';
project = 'Dl-Ven_snaBAC-mCh_v4';
DataPath = [DropboxFolder 'ProcessedEnrichmentData/' project '/'];
% Params
K = 3;
w = 7;
FigPath = ['C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\Dl-Ven_snaBAC-mCh_v4\input_output_traces\'];
mkdir(FigPath);
% intermediate HMM results
load([DataPath 'nucleus_struct_protein.mat'])    


%% specify plotting constraints
PixelSize = nucleus_struct_protein(1).PixelSize;
max_jump = 0.8 / 21; %um
max_nc_dist = 2.5; %um
min_n = 75;
n_vec = [nucleus_struct_protein.N];
n_ids = find(n_vec>=min_n);

keep_ids = [];
for n = n_ids
    xn_vec = nucleus_struct_protein(n).xPos;
    yn_vec = nucleus_struct_protein(n).yPos;
    xp_vec = nucleus_struct_protein(n).xPosParticle - xn_vec;
    yp_vec = nucleus_struct_protein(n).yPosParticle - yn_vec;
    
    t_vec = nucleus_struct_protein(n).time;
    dt_vec = diff(t_vec);
    
    nc_dist_vec = sqrt(xp_vec.^2 + yp_vec.^2)*PixelSize;
    
    r_vec = sqrt(diff(xp_vec).^2 + diff(yp_vec).^2)*PixelSize;
    
    if all(r_vec./dt_vec<=max_jump) && all(nc_dist_vec<=max_nc_dist)
        keep_ids = [keep_ids n];
    end
end


%%

rng(341);
cmap = brewermap(9,'Set2');

for p = keep_ids
    
   
    % extract basic time, fluo, and protein fields
    time = nucleus_struct_protein(p).time/60;    
    
    fluo = nucleus_struct_protein(p).fluo;    
    
    spot_protein = nucleus_struct_protein(p).spot_protein_vec;    
    serial_protein = nucleus_struct_protein(p).serial_null_protein_vec;    



    % make figure
    burst_surge_fig = figure('Visible','off');
    hold on       
    yyaxis right
%     p_pt = area(time,spot_protein,'FaceColor',cmap(2,:),'EdgeAlpha',0,'FaceAlpha',0.25);
    plot(time,spot_protein,'Color',cmap(2,:),'LineWidth',1.5);
    plot(time,serial_protein,'-','Color',cmap(3,:),'LineWidth',1.5);
    ax = gca;
    ax.YColor = cmap(2,:);
    ylabel('Dorsal al locus (detrended)')
    
    yyaxis left
    % fluo trend
    area(time,fluo,'FaceColor',[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0);%,'LineWidth',1);

    % plot inferred burst events
    ax = gca;
    xlim([nanmin(time) nanmax(time)])
    y_lim = ax.YLim;
    ylim([0 max(y_lim)])
    y_lim = ax.YLim;   
    ax.YColor = 'black';
        
    ylabel('MS2 spot intensity')
    
    legend('spot fluorescence','protein at locus','protein at control')    
    xlabel('minutes into nc14')
    set(gca,'Fontsize',14);
    box on;
    
    saveas(burst_surge_fig,[FigPath 'trace_pt' num2str(p) '.png'])
    saveas(burst_surge_fig,[FigPath 'trace_pt' num2str(p) '.pdf'])
    close all
        
end