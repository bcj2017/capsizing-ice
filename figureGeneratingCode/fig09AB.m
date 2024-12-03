clear; close all;
expname = '2022-07-13-a';
t = figure('Units','Normalized','Position',[0 0 .6 .7]);

%% ------------------------------------------------------------------------------------
%% Plot settings
%% ------------------------------------------------------------------------------------
lw = 2.5;
cmap = cool;
clr_idx = 115;
shapeshade_idx = 40; shapeshade_int = 60;
fontname = 'times'; fontsize = 25;
kk = 1e-2; ticklengthfac = 2.5; ticklinewidthfac = 150;
upfac = .05;

%% ------------------------------------------------------------------------------------
%% Chosen frames, path to data set
%% ------------------------------------------------------------------------------------
ts_arr = [78, 346, 746]; 
pathtodata = ['../data/',expname,'/'];

%% ------------------------------------------------------------------------------------
%% Plot for 3 chosen frames
%% ------------------------------------------------------------------------------------
curv_shift = 2750; 
for j = 1:length(ts_arr) % iterate over chosen frames
    % Load data
    ts = ts_arr(j); %set time
    load([pathtodata,num2str(ts),'.mat']);
    % Shift curvature function to display peaks more clearly
    curvature = [curvature(curv_shift:end); curvature(1:curv_shift)] / mean(curvature);

    %% Plot 1: Curvature versus arclength
    subplot(3,2,2*(j-1)+1)
    s = linspace(0,1,length(curvature));
    pl = plot(s,curvature); pl.LineWidth = lw; pl.Color = cmap(clr_idx,:);
    if j == 2 % y axis label
        hYLabel = ylabel('curvature, $\kappa/\bar{\kappa}$','Interpreter','latex', 'Position', [-0.1, 0.5, 0]); 
        set(hYLabel, 'FontSize', fontsize, 'VerticalAlignment', 'middle');
        hYLabel.Position(1) = hYLabel.Position(1) - 0.05;
        hYLabel.Position(2) = hYLabel.Position(2) + 2; 
    end
    % Set xticks for top two plots
    xticks([0:.2:1]); xticklabels([]);
    if j == 3 % and for bottom plot
        xlabel('arclength, $s/L$','Interpreter','latex','FontName',fontname,'FontSize',fontsize)
        xticklabels([0 0.2 0.4 0.6 0.8 1]);
    end
    % Plot details
    xlim([0,1]); ylim([0,5])
    ax = gca; ax.FontName = fontname; ax.FontSize = fontsize; 
    ax.TickLength = ticklengthfac*kk*ones(1,2); ax.LineWidth = ticklinewidthfac*kk; % Make tick marks longer and thicker.
    % Move up
    pos = ax.Position; pos(2) = pos(2) + upfac * (j-1); ax.Position = pos;

    %% Plot 2: Fourier transform
    subplot(3,2,2*(j-1)+2); hold on; box on
    ylimit = 21; % maximum value we plot over
    % Compute fourier transform and |a_k|^2/N
    ft = abs(fft(curvature)); ft = ft.^2 / length(x_interp)^2;
    pl = plot(1:1:ylimit-1,ft(2:ylimit));  
    xticks([0:5:ylimit]);
    pl.LineWidth = lw; pl.Color = cmap(clr_idx,:);
    if j == 2
        hYLabel2 = ylabel('Fourier amplitude, $|a_k|^2/N^2$','Interpreter','latex');%, 'Position', [-0.1, 0.5, 0],'FontSize',fontsize);
        % set(hYLabel2, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
        % hYLabel2.Position(1) = hYLabel2.Position(1) - 3.5;
        % hYLabel2.Position(2) = hYLabel2.Position(2) + 155; 
    end
    % Set xtick labels to only be on bottom plot:
    xticklabels([]);
    if j == 3
        xlabel('$k$','Interpreter','Latex','FontName',fontname,'FontSize',fontsize)
        xticklabels([0; 5; 10; 15]);
    end

    ax = gca; ax.FontName = fontname; ax.FontSize = fontsize; 
    ax.TickLength = ticklengthfac*kk*ones(1,2); ax.LineWidth = ticklinewidthfac*kk; % Make tick marks longer and thicker.
    pos = ax.Position; pos(1) = pos(1)+.02; pos(2) = pos(2) + upfac * (j-1); ax.Position = pos;
    xlim([0,16]); ylim([0,.1]); 

    %% ------------------------------------------------------------------------------------
    %% Plot shape in corner of Fourier
    %% ------------------------------------------------------------------------------------
    subplotPos = ax.Position;
    % Set up subplots for each polygon, plotting shapes to relative scale.
    if j == 1
        axes('Position', [subplotPos(1) + 0.44 * subplotPos(3), subplotPos(2) + 0.32 * subplotPos(4), .52*subplotPos(3), .62*subplotPos(4)]);
    elseif j == 2
        axes('Position', [subplotPos(1) + 0.47 * subplotPos(3), subplotPos(2) + 0.28 * subplotPos(4), .52*subplotPos(3), .62*subplotPos(4)]);
    else
        axes('Position', [subplotPos(1) + 0.5 * subplotPos(3), subplotPos(2) + 0.28 * subplotPos(4), .52*subplotPos(3), .62*subplotPos(4)]);
    end

    % Plot polyshape and boundary
    hold on
    ss = polyshape(x_interp,y_interp);
    vertices = ss.Vertices;
    fill(vertices(:,1), vertices(:,2), cmap(shapeshade_idx,:), 'FaceAlpha', .6); %1 is opaque, 0 is invisible
    pl = plot(x_interp,y_interp); pl.LineWidth = lw; pl.Color = 'k';
    if j == 1 % set x limits for consistency / scale
        ymin = min(y_interp); ymax = max(y_interp);
    end
    ylim([ymin,ymax]);
    axis off; axis equal
    shapeshade_idx = shapeshade_idx + shapeshade_int;
end

%% ------------------------------------------------------------------------
%% Save file
%% ------------------------------------------------------------------------
filename = 'fig09ABdraft.pdf';
exportgraphics(t,filename,'ContentType','vector')
