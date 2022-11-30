clear all
close all
clc

%% Part 1: Read data
load('./data/2020-02-23-optoknirps_new_embryo6_lin.mat');
load('./data/CompiledParticles.mat')

% frames to be processed
start_frame = 40;
final_frame = 149;

%% Part 2: Quality check on nuclei tracking

% Quality check and save the nuclei that passed the test
% Requires the nuclei to have continuous trace until lateral movement
X_pass = [];
Y_pass = [];

X_fail = [];
Y_fail = [];

sch_pass = [];

for i = 1:size(schnitzcells,2)
    index_s = find(schnitzcells(i).frames == start_frame);
    index_f = find(schnitzcells(i).frames == final_frame);
    fluo_temp = max(schnitzcells(i).Fluo(index_s:index_f,:),[],2);
    result = sum(isnan(fluo_temp));
    % quality check
    if ~isempty(index_s) && ~isempty(index_f) && (index_f-index_s == final_frame-start_frame)  ...
            && (result == 0)
        sch_pass = [sch_pass, i];
        X_pass = [X_pass, schnitzcells(i).cenx(index_s)];
        Y_pass = [Y_pass, schnitzcells(i).ceny(index_s)];
    else
        X_fail = [X_fail, schnitzcells(i).cenx(index_s)];
        Y_fail = [Y_fail, schnitzcells(i).ceny(index_s)];
    end
end

figure(1)
plot(X_pass,Y_pass,'o')
hold on
plot(X_fail,Y_fail,'o')
axis equal
xlim([0 1024])
ylim([0 256])
title('Nuclei that passed quality check')


%% Part 3: Compile the schnitzcells together

% Initialize storage
processed_data(1).xcoord = [];
processed_data(1).ycoord = [];

processed_data(1).schnitznum = [];
processed_data(1).NuclearFluo = [];
processed_data(1).SpotFluo = [];

% Compile all the nuclei in each frame and assign basic info
for i = 1:length(sch_pass)
    sch_num = sch_pass(i);
    index_s = find(schnitzcells(sch_num).frames == start_frame);
    index_f = find(schnitzcells(sch_num).frames == final_frame);
    
    for j = index_s:index_f
        sch_now = sch_pass(i);
        frame_now = schnitzcells(sch_now).frames(j);
        x_coord = schnitzcells(sch_now).cenx(j);
        y_coord = schnitzcells(sch_now).ceny(j);
        fluo = max(schnitzcells(sch_now).Fluo(j,:));
        try
            processed_data(frame_now).xcoord = [processed_data(frame_now).xcoord, x_coord];
            processed_data(frame_now).ycoord = [processed_data(frame_now).ycoord, y_coord];
            processed_data(frame_now).schnitznum = [processed_data(frame_now).schnitznum, sch_now];
            processed_data(frame_now).SpotFluo = [processed_data(frame_now).SpotFluo, 0];
            processed_data(frame_now).NuclearFluo = [processed_data(frame_now).NuclearFluo fluo];  
        catch
            processed_data(frame_now).xcoord = x_coord;
            processed_data(frame_now).ycoord = y_coord;
            processed_data(frame_now).schnitznum = sch_now;
            processed_data(frame_now).SpotFluo = 0;
            processed_data(frame_now).NuclearFluo = fluo;
        end
    end
end

% Assign particle info to all the nuclei
for i = 1:size(CompiledParticles{1,1},2)
    schnitz_num = CompiledParticles{1,1}(i).schnitz;
    for j = 1:size(CompiledParticles{1,1}(i).Frame,2)
        frame = CompiledParticles{1,1}(i).Frame(j);
        if frame<=final_frame
            num = find(processed_data(frame).schnitznum==schnitz_num);
            processed_data(frame).SpotFluo(num) = CompiledParticles{1,1}(i).Fluo(j);
        end
    end
end

%{
frame_plot = 149;
xpos = processed_data(frame_plot).xcoord;
ypos = processed_data(frame_plot).ycoord;

% try to plot and test a bit...
pts = [xpos' ypos'];

fig = figure(2);
[v,c] = voronoin(double(pts));

for i = 1:length(c)
    if all(c{i}~=1)
    x = v(c{i},1);
    y = v(c{i},2);
    %a = processed_data(frame_plot).SpotFluo(i);
    a = processed_data(frame_plot).NuclearFluo(i);
    patch(x,y,a);
    colorbar
    caxis([0 9E5])
    end
end

axis equal
xlim([0 1024])
ylim([0 256])
%}

% Convert the format: Store them into individual nuclei traces
nuclei_fluo_traces = zeros(length(sch_pass),final_frame);
spot_fluo_traces = zeros(length(sch_pass),final_frame);
x_pos = zeros(length(sch_pass),final_frame);
y_pos = zeros(length(sch_pass),final_frame);

for i = 1:length(sch_pass)
    for j = start_frame:final_frame
        nuclei_fluo_traces(i,j) = processed_data(j).NuclearFluo(i);
        spot_fluo_traces(i,j) = processed_data(j).SpotFluo(i);
        x_pos(i,j) = processed_data(j).xcoord(i);
        y_pos(i,j) = processed_data(j).ycoord(i);
    end
end

%% Part 4: Visualize input and output

factor = 1E4;

%for i = 140
for i = start_frame:final_frame
    % Plot voronoi first
    xpos = processed_data(i).xcoord;
    ypos = processed_data(i).ycoord;

    pts = [xpos' ypos'];

    fig = figure(2);
    clf
    [v,c] = voronoin(double(pts));

    for j = 1:length(c)
        if all(c{j}~=1)
        x = v(c{j},1);
        y = v(c{j},2);
        a = processed_data(i).NuclearFluo(j); % plot nuclear fluorescence
        patch(x,y,a);
        colorbar
        caxis([1E5 8E5])
        colormap('jet');
        end
    end

    axis equal
    xlim([0 1024])
    ylim([0 256])
    hold on
    
    % Plot Spots
    for j = 1:length(processed_data(i).xcoord)
        x = processed_data(i).xcoord(j);
        y = processed_data(i).ycoord(j);
        marker_size = processed_data(i).SpotFluo(j)/factor;
        if marker_size>0
            plot(x, y, '.', 'MarkerSize',marker_size,'Color','#A2142F')
        end
    end
    hold off
    title(['Time: ',num2str(ElapsedTime(i)),' min'])
    
    saveas(fig,['./figure/test/Frame_',num2str(i),'.jpg'])
    
end
