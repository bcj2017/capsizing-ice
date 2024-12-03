%% Generate figure 5 of data
clear; close all
addpath('functions'); % Add functions directory to path

%% Theoretical parameters
theoryalpha = 79.85 * pi / 180; % Convert to radians
theoryN = 2 * pi / theoryalpha; % Theoretical value for N
numrows = 2; numcols = 5; % Number of rows and columns in subplot

%% Colormap settings
cmap = cool; % Cool colormap
startclridx = 20; clridx = startclridx; % Start color index
clrint = 22; % color interval
simcolor = [0.2 0.2 0.6]; % Color for simulation data
cpgray = gray; gray_idx = 170; % gray coloring
% Arrow settings:
arrow_lw = 2.5; arrowcolor = cpgray(gray_idx,:); 
eps = .3; tri_size = 0.7; y_arrowspacing = .7; 
theta_start = pi/3-2*eps; theta_end = pi/2-eps;
theta = linspace(theta_start, theta_end, 100);

%% Generic plot settings
t = figure('Units', 'normalized', 'Position', [0 0 1 1]); % Create figure
lredfac = 1.5; % Factor to reduce marker size in insets
simmarker = 'o'; expmarker = '.';
expscatsize = 80; simmarkersize = 8; % Marker size
lw = 1.8; trans = 0.8; fontname = 'times'; fontsize = 15; % Other plot parameters
% Tick length + width parameters:
k = 1e-2; ticklengthfac = 2.5; ticklinewidthfac = 150;
lwalpha = 2; % linewidth for alpha dotted line
desired_aspect_ratio = 1/3; % for histogram plot window sizing
% Parameters for spacing subplots:
xspacing = 8e-3; histspacing = 2.5e-2;
maxflip = 12; % Last flip included for simulation data

%% Path to data and list of experiments
pathtoexp = '../data/experimental-measured-quantities/'; % Path to experimental data
pathtosim = '../data/simulationData/'; % Path to simulation data
S = dir(fullfile(pathtoexp, '**', '*.mat')); names = {S.name}; % File names
totalSimTime = 27.97;

%% ------------------------------------------------------------------------
%% Panel A1: area plot
%% ------------------------------------------------------------------------
h1 = subplot(numrows, numcols, 1); hold on;

% Plot power law line
tt = linspace(0, 1, 1000);
plot(tt, (1 - tt).^(8 / 5), '--', 'Color', 'k', 'Linewidth', lw);

% Plot experimental data for area
for i = 1:length(names) % iterate through experiments we have data on
    % Load data for experiment
    matname = char(names(i));
    load([pathtoexp,matname]);
    
    areaTimes = Area(:,1); 
    Area = Area(:,2) / Area(1,2); % normalize area by initial data point
        
    % Polyfit to find final time
    xx = polyfit(areaTimes,Area.^(5/8),1);
    tf = -xx(2)/xx(1); 
    areaTimes = areaTimes/tf;
    
    % Plot
    p = scatter(areaTimes,Area,expscatsize,cmap(clridx,:),'filled', ...
        'MarkerFaceAlpha', trans, 'MarkerEdgeAlpha', 1, ...
        'LineWidth',lw);

    % Increase color index to distinguish between experiments
    clridx = clridx+clrint;
end

% Plot simulation data for area
load([pathtosim,'simulationData.mat']);
sim_times = Area(:,1) / totalSimTime; areas = Area(:,2);
p = plot(sim_times,areas,simmarker); 
p.LineWidth = lw; p.MarkerSize = simmarkersize; p.Color = simcolor;

% Plot details
xlabel('time, $t/t_f$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize); 
ylabel('area, $A/A_0$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize); 
xlim([0,1]); ylim([0,1]); 
box on; axis equal; axis square; 
ax = gca; 
ax.FontName = fontname; ax.FontSize = fontsize; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = ticklinewidthfac*k; 

clridx = startclridx; % reset color index for next plot

%% ------------------------------------------------------------------------
%% Panel A2: log area plot
%% ------------------------------------------------------------------------
% Define inset positioning for log-log plot:
subplotPos = get(h1, 'Position');
axes('Position', [subplotPos(1) + 0.53 * subplotPos(3), subplotPos(2) + 0.45 * subplotPos(4), .41*subplotPos(3), .41*subplotPos(4)]);

% Plot power-law
epsilon = 1e-2; tt = linspace(0,1-epsilon,1000);
loglog(1-tt,(1-tt).^(8/5), '--', 'Color','k','Linewidth',lwalpha);
hold on;

for j = 1:length(names) %iterate through experiments we have data on
    % Load data for experiment
    matname = char(names(j));
    load([pathtoexp,matname]);
    
    areaTimes = Area(:,1); 
    Area = Area(:,2) / Area(1,2); % normalize area by initial data point
        
    % Polyfit to find final time
    xx = polyfit(areaTimes,Area.^(5/8),1);
    tf = -xx(2)/xx(1);
    areaTimes = areaTimes/tf;
        
    % Plot
    p = scatter(1-areaTimes,Area,expscatsize/lredfac,cmap(clridx,:),'filled', ...
        'MarkerFaceAlpha', trans, 'MarkerEdgeAlpha', 1, ...
        'LineWidth',lw);

    % Increase color index to distinguish between experiments
    clridx = clridx+clrint; 
end

% Simulation plot
p = plot(1-sim_times,areas,simmarker); 
p.LineWidth = lw; p.MarkerSize = simmarkersize; p.Color = simcolor;

% Plot details
xlabel('$1-t/t_f$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize); 
ylabel('$A/A_0$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize); 
box on; axis equal; axis square; ax = gca; 
ax.FontSize = fontsize; ax.FontName = fontname; ax.TickLength = 2*ticklengthfac*k*ones(1,2); ax.LineWidth = 1*ticklinewidthfac*k; 
% Define specific tick locations
xTickLocations = logspace(-1, 2, 5); yTickLocations = logspace(-2, 4, 6); 

clridx = startclridx; % reset color index for next plot

%% ---------------------------------------------------------------------------
%% Panel B: shape mode plot
%% ---------------------------------------------------------------------------
plotShapeMode = subplot(numrows,numcols,2); hold on
pos = plotShapeMode.Position; pos(1) = pos(1) - .008; plotShapeMode.Position = pos; % adjust position

for i = 1:length(names) %iterate through experiments we have data on
    matname = char(names(i));
    load([pathtoexp,matname]);

    % after-flip data only:
    shapeModeTimes = shapeMode(2:end,1); shapeModeTimes = shapeModeTimes / 60;
    shapeMode = shapeMode(2:end,2);

    p = scatter(shapeModeTimes,shapeMode,expscatsize,cmap(clridx,:),'filled', ...
        'MarkerFaceAlpha', trans, 'MarkerEdgeAlpha', 1, ...
        'LineWidth',lw);

    clridx = clridx+clrint; % adjust color for next experiment
end
%% Plot simulation shapemode data
load([pathtosim,'simulationData.mat']);
sim_times = shapeMode(2:end,1); sim_modes = shapeMode(2:end,2); 
p = plot(sim_times,sim_modes,simmarker);
p.LineWidth = lw; p.MarkerSize = simmarkersize; p.Color = simcolor;

% Theoretical value:
x = linspace(0,30,100);y = theoryN*ones(1,100);
plot(x,y,'--k','Linewidth',lw);

% Plot details:
xlabel('time, $t$ (min)','Interpreter','latex'); ylabel('shape mode, $k_*$','Interpreter','Latex');
xlim([0,30]); ylim([0,15]);
box on; axis square; ax = gca; 
ax.FontName = fontname; ax.FontSize = fontsize; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = ticklinewidthfac*k; 

%% Set up right y-axis
yyaxis right;
ylim([0,15]);
yticks(0:5:15); % Set the same tick positions
custom_labels_right = {'','$2\pi/5$','$\pi/5$','$2\pi/15$'};
yticklabels(custom_labels_right);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.YColor = 'k'; % Change the color of the y-axis line and ticks to black

clridx = startclridx; % reset color index for next plot

%% ------------------------------------------------------------------------
%% Panel B2: curvature inset
%% ------------------------------------------------------------------------
% Inset positioning:
subplotPos = get(plotShapeMode, 'Position');
axes('Position', [subplotPos(1) + 0.6 * subplotPos(3), subplotPos(2) + 0.485 * subplotPos(4), .35*subplotPos(3), .4*subplotPos(4)]);

%% Inset: load and plot data for chosen frame
load('sample-data/fig5curvature.mat'); % load
p = plot(s,ccva); % plot curvature versus arclength
xlim([min(s),max(s)]); ylim([-1,3]); % plot settings
% Plot details:
p.LineWidth = lw; p.Color = cmap(115,:);
box on; axis square; ax = gca; 
ax.FontName = fontname; ax.FontSize = fontsize; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = ticklinewidthfac*k; 
% Custom xtick labels
ax.XTick = [min(s),max(s)];
customLabels = {'0', '1'}; ax.XTickLabel = customLabels;
ax.TickLabelInterpreter = 'latex';

% Set x and y labels and reposition to fit inset better
h = ylabel('$\kappa/\bar{\kappa}$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize);
pos = get(h, 'Position'); pos(1) = pos(1) + 1.6; set(h, 'Position', pos); 
h = xlabel('$s$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize);
pos = get(h, 'Position'); pos(2) = pos(2) + .7; set(h, 'Position', pos); 

% Mark peaks on curvature plot
plotLineWithVerticalMarkers(pt1, pt2, height, lw); 

%% ------------------------------------------------------------------------
%% Panel B2: shape boundary inset
%% ------------------------------------------------------------------------
% Set new axis positions:
axes('Position', [subplotPos(1) + 0.23 * subplotPos(3), subplotPos(2) + 0.59 * subplotPos(4), .26*subplotPos(3), .26*subplotPos(4)]);
hold on

% Plot boundary:
ss = polyshape(ex_shift,ey_shift); vertices = ss.Vertices;
fill(vertices(:,1), vertices(:,2), cmap(115,:), 'FaceAlpha', .7);
p = plot(ex_shift,ey_shift); p.LineWidth = lw; p.Color = 'k';

%% Plot arrowhead moving counterclockwise:
r_arrow = mean(sqrt((ex_shift).^2 + (ey_shift).^2))*1.2; % Pick radius for arrow
y_arrowheight = max(ey_shift) + y_arrowspacing - r_arrow;
% Compute x and y-coordinates of arrow body:
x_arrow = r_arrow * cos(theta);
y_arrow = r_arrow * sin(theta) + y_arrowheight;
% Plot arrow body:
plot(x_arrow,y_arrow,'Linewidth',arrow_lw,'Color',arrowcolor)
% Details for arrowhead:
center_x = x_arrow(end); center_y = y_arrow(end); 
theta_arrowhead = theta_end; 
% Plot arrowhead:
drawCurvedArrow(theta_end, theta_start, tri_size, center_x, center_y, theta_arrowhead, arrowcolor);

axis square; axis off; axis equal
%% ---------------------------------------------------------------------------
%% Panel C2: histogram of flip angles
%% ---------------------------------------------------------------------------

allRotationAngles = []; % initialize arrays to store all data for experimental histogram:
for i = 1:length(names) %iterate through experiments
    % Load data
    matname = char(names(i));
    load([pathtoexp,matname]);

    % Append data into array for histogram
    allRotationAngles = [allRotationAngles; rotationAngle(:,2)];
end

%% Flip angle vs. time: subplot positioning
ax = subplot(numrows,numcols,3);
pos1 = ax.Position; 

%% Histogram plot
% Subplot and positioning + sizing
ax = subplot(numrows,numcols,4); hold on
pos2 = ax.Position;
pos2(1) = pos2(1) + xspacing - histspacing;  pos2(2) = .637;  pos2(3) = pos2(3) * desired_aspect_ratio; pos2(4) = pos2(4) * .687;
ax.Position = pos2;

% Compute histogram data
nbins = 40;
[counts, edges] = histcounts(allRotationAngles,nbins);
% Calculate bin centers
bin_centers = edges(1:end-1) + diff(edges)/2;
% Create horizontal bar plot
p = barh(bin_centers, counts);
p.BarWidth = 1; p.FaceColor = cmap(115,:); p.EdgeColor = 'none'; 
xlim([0,20]);

% Plot integer multiples of theoretical alpha:
x = linspace(0,20,100);y = 2*theoryalpha*ones(1,100);
for j = 1:5
    if y(1) ~= 0
        plot(x,y,'--k','Linewidth',lw);
    end
    y = y - theoryalpha;
end
uistack(p,'top'); % bring histogram to top
% Plot details:
xlabel('counts','Interpreter','latex'); box on
ylim([-pi pi]); yticks([-pi -pi/2 0 pi/2 pi]);
ax = gca; 
ax.FontName = fontname; ax.FontSize = fontsize; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = ticklinewidthfac*k; 
ax.YTickLabel = [];
%% ------------------------------------------------------------------------
%% Panel C1: flip angle versus time plot
%% ------------------------------------------------------------------------
subplot(numrows,numcols,3); hold on
% Plot integer multiples of theoretical alpha:
x = linspace(0,30,100); y = 2*theoryalpha*ones(1,100);
for j = 1:5
    if y(1) ~= 0
        plot(x,y,'--k','Linewidth',lw);
    end
    y = y - theoryalpha;
end

% Plot experimental flip angle data
for i = 1:length(names) % iterate through experiments
    % Load data for each experiment
    matname = char(names(i));
    load([pathtoexp,matname]);   
    
    flipTimes = flipTimes(1:length(rotationAngle)) ./ 60; 
    p = scatter(flipTimes,rotationAngle,expscatsize,cmap(clridx,:),'filled', ...
            'MarkerFaceAlpha', trans, 'MarkerEdgeAlpha', 1);
    clridx = clridx + clrint; % increase color index to distinguish between experiments
end
% Plot simulation data
load('../data/simulationData/simulationData.mat');
sim_fliptimes = rotationAngle(1:maxflip,1); sim_flipangle = rotationAngle(1:maxflip,2);
p = plot(sim_fliptimes,sim_flipangle,simmarker);
p.LineWidth = lw; p.MarkerSize = simmarkersize*1.2; p.Color = simcolor;

% Plot details:
xlim([0,30]);
xlabel('time, $t$ (min)','Interpreter','latex'); ylabel('flip angle, $\Delta \phi$','Interpreter','Latex')
box on; axis square; ax = gca; 
ylim([-pi pi]); yticks([-pi -pi/2 0 pi/2 pi]); yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yyaxis right; ax2 = gca(); 
ax.YAxis(1).Color = ax.XAxis.Color; ax.YAxis(2).Color = ax.XAxis.Color; ax.YAxis(2).TickValues = []; 
ax.FontName = fontname; ax.FontSize = fontsize; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = ticklinewidthfac*k; 
pos1(1) = pos1(1) + xspacing; ax.Position = pos1;

clridx = startclridx; % restart color index

%% ---------------------------------------------------------------------------
%% Panel D1: flip intervals versus time plot
%% ---------------------------------------------------------------------------
% Set up subplot and position
h = subplot(numrows,numcols,8); hold on
posfliptime = h.Position; posfliptime(2) = posfliptime(2)+.18; h.Position = posfliptime;

deltaFlipTimesAll = []; % initialize array of fliptimes for histogram
for i = 1:length(names) %iterate through experiments
    % Load experimental data
    matname = char(names(i));
    load([pathtoexp,matname]);
    
    % Reorganize flip interval data:
    flipTimes = flipTimes / 60; % convert to minutes
    deltaFlipTimes = flipTimes(2:end) - flipTimes(1:end-1);

    % Save all flip intervals for histogram
    for j = 1:length(deltaFlipTimes)
        deltaFlipTimesAll(end+1) = deltaFlipTimes(j);
    end
    
    % Scatter plot flip intervals versus time
    p = scatter(flipTimes(1:end-1),deltaFlipTimes,expscatsize,cmap(clridx,:),'filled', ...
        'MarkerFaceAlpha', trans, 'MarkerEdgeAlpha', 1, ...
        'LineWidth',lw);

    clridx = clridx+clrint; % increase color index to distinguish between experiments
end
% Plot simulation flip intervals
delta_sim_fliptimes = sim_fliptimes(2:end) - sim_fliptimes(1:end-1);
p = plot(sim_fliptimes(1:end-1),delta_sim_fliptimes,'o');
p.LineWidth = lw; p.MarkerSize = simmarkersize; p.Color = simcolor;

% Plot details:
ylim([0,6]);
xlabel('time, $t$ (min)','Interpreter','latex'); ylabel('flip interval, $\Delta t$ (min)','Interpreter','latex')
box on; axis square; ax = gca; 
ax.FontSize = fontsize; ax.FontName = fontname; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = 1*ticklinewidthfac*k; % Make tick marks thicker.
posD1 = ax.Position; posD1(1) = posD1(1) + xspacing; ax.Position = posD1; % reposition

%% ---------------------------------------------------------------------------
%% Panel D1: histogram of flip intervals
%% ---------------------------------------------------------------------------
% Set up subplot and positioning
ax = subplot(numrows,numcols,9);
pos2 = get(gca,'Position'); % [left bottom width height]
pos2(1) = pos2(1) + xspacing - histspacing; pos2(2) = pos2(2) + .2333; pos2(3) = pos2(3) * desired_aspect_ratio; pos2(4) = pos2(4) * .687;
ax.Position = pos2;

% Compute histogram data
nbins = 20;
[counts, edges] = histcounts(deltaFlipTimesAll,nbins);
% Calculate bin centers
bin_centers = edges(1:end-1) + diff(edges)/2;
% Create horizontal bar plot
p = barh(bin_centers, counts); % Customize color as needed
set(p, 'BarWidth', 1); % Set bar width to 1 to fill the space
p.FaceColor = cmap(115,:); p.EdgeColor = 'none'; 

% Plot details:
xlim([0,20]); ylim([0,6]);
box on
ax.FontSize = fontsize; ax.FontName = fontname; ax.TickLength = ticklengthfac*k*ones(1,2); ax.LineWidth = 1*ticklinewidthfac*k; % Make tick marks thicker.
ax.YTickLabel = [];
xlabel('counts','Interpreter','latex')

%% ------------------------------------------------------------------------
%% Save file
%% ------------------------------------------------------------------------
filename = 'fig05ABEFdraft.pdf';
% exportgraphics(t,filename,'ContentType','vector')
