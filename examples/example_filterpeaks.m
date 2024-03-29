%% example for filtering peaks from signal traces

% load some signal trace data
load('./examples/example_ROItraces.mat')


%% visualize traces
figure
PlotAxes =gca;
plotTraces(ROIs,PlotAxes)

%% Pick out peaks using some cutoffs
options = struct();

% polynomial fit order for detrending
options.detrendorder = 10;

% svg order and window size for smoothing signal to get curvature
options.order = 3;
options.framelen = 15;

% scaling factor for universal threshold
% for defining peak in raw data
options.cutscalepeak = 0.8;
% for defining sufficient negative curvature in smoothed data
options.cutscalecurvature = 0.5;

[havepeakind,ROIs] = filterForPeaks(ROIs,options)

% plot resulting peaks
plotTraces(ROIs,PlotAxes)
