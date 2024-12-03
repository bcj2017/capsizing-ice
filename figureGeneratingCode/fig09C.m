clear; close all;
addpath('functions')
t = figure('Units','Normalized','Position',[0 0 .6 1/3]);
expname = '2022-07-13-a';
exptime = 29*60;

%% ------------------------------------------------------------------------------------
%% Plot settings
%% ------------------------------------------------------------------------------------
lw = 3; fontname = 'times'; fontsize = 25; 

%% ------------------------------------------------------------------------------------
%% Plot kymograph
%% ------------------------------------------------------------------------------------
% Generate kymograph:
[h,cb] = generateKymograph(expname,exptime);

% Color bar settings and title
cb_pos = cb.Position; cb.FontName = fontname; cbAxes = cb.Axes;
clim([0,.1]); % axis limits

cb.Label.String = 'Fourier amplitude, $|a_k|^2/N^2$'; cb.Label.Interpreter = 'latex';

%% ------------------------------------------------------------------------------------
%% Plot shape modes versus time in white line
%% ------------------------------------------------------------------------------------
% Load and plot shape mode data
load(['../data/experimental-measured-quantities/',expname,'.mat']);
% White line shape mode plot:
p = plot3(shapeMode(2:end,2),shapeMode(2:end,1)/60,zeros(size(shapeMode(2:end,1))));
p.Color = 'w'; p.LineWidth = lw; % line settings
% Scatter plot on top:
scatter3(shapeMode(2:end,2),shapeMode(2:end,1)/60,zeros(size(shapeMode(2:end,1))),70,'w','filled');

% Plot settings
ax = gca; ax.FontName = fontname; ax.FontSize = fontsize; 
pos = ax.Position; pos(3) = .75*pos(3); ax.Position = pos; % reposition

%% ------------------------------------------------------------------------
%% Save file
%% ------------------------------------------------------------------------
filename = 'fig09Cdraft.pdf';
% exportgraphics(t,filename,'ContentType','vector')
