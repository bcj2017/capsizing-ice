function [h,cb] = generateKymograph(expname,exptime)

%% ------------------------------------------------------------------------
%% Import data
%% ------------------------------------------------------------------------
pathtobdies = ['../data/',expname,'/'];
pathtodata = ['../data/experimental-measured-quantities/',expname,'.mat'];
load(pathtodata);
ts_arr = shapeMode(:,1);

%% ------------------------------------------------------------------------
%% Compute Fourier amplitudes for each time stamp
%% ------------------------------------------------------------------------
numframes = length(ts_arr); % number of iterations
fourierAmplitudes = cell(1,numframes); % cell containing all Fourier amplitudes at each iteration
for j = 1:numframes % iterate over frames to compute Fourier amplitudes
    ts = ts_arr(j);
    % Load boundary for this time stamp:
    load([pathtobdies,num2str(ts),'.mat']);
    n = length(x_interp);
    % Normalize and mean-zero curvature plot:
    curvature = curvature / mean(curvature); % normalize for kymograph 
    curvature = curvature - mean(curvature); % mean zero (to ignore k=0 mode)

    % Fourier transform of curvature, extract |a_k|^2 / n:
    ft = abs(fft(curvature)).^2 / n;
    % Take number of modes we plot:
    nummodes = 18;
    ft = ft(1:nummodes); % for kymograph

    % Save data in cell fourier_amplitudes
    fourierAmplitudes{j} = ft;
end % end iteration over frames

%% ------------------------------------------------------------------------
%% Reorder fourier_amplitudes cell to plot kymograph
%% ------------------------------------------------------------------------
shapeModeAmplitudes = zeros(nummodes-1, numframes+1); % delete 0 mode --> nummodes-1
for i = 2:nummodes 
    for j = 1:numframes
        shapeModeAmplitudes(i-1,j) = fourierAmplitudes{j}(i);
    end
    shapeModeAmplitudes(i-1,end) = fourierAmplitudes{end}(i);
end
% Add extra dummy point to extend to end of experiment
ts_arr(end+1) = exptime;

%% ------------------------------------------------------------------------
%% Plot kymograph
%% ------------------------------------------------------------------------
k = .5:1:nummodes-1.5; % to center k's in middle of actual value
tmin = ts_arr / 60; % convert to minutes
[X,Y] = meshgrid(k,tmin);
h = surf(X,Y,shapeModeAmplitudes'); % Plot kymograph
hold on
ax = gca;
% Reposition:
position1 = ax.Position;
position1(3) = 1.2*position1(4); position1(4) = position1(4)*.9; position1(2) = position1(2)+.06;
ax.Position = position1;
view(90, -90); % bird's-eye view
xlim([min(k),max(k)]); ylim([min(tmin),max(tmin)]); % x and y limits
xticks([5,10,15,20]);
xlabel('$k$','Interpreter','latex'); ylabel('time (min)','Interpreter','latex')
colormap("cool");
cb = colorbar;

% Turn off the mesh lines
set(h, 'EdgeColor', 'k', 'LineWidth', .15); % 'k' for black, adjust LineWidth as needed
end