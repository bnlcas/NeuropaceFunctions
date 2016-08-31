function [] = plot_spectrogram_rns_evnt(evnt, appelation, smoothing, varargin)
%% Takes an evnt structure and loads the 4 channels of RAW data for to generate spectrograms
% For each Channel:
% 
% Inputs:
%
% evnt - event structure - needed to load directories for the raw data
%
% appelation - title to be used on figures;
%
% smoothing - smoothing parameter on data
%
% varargin:
%
% input 1: y axis range setting
%
% input 2: normalize plots (default set to normalize to power in the 0 -
% 125 hz range
% input 3: hold plots - parameter to determine whether to plot new plots of
% overwrite plots.

scale_y = false;
normalize = false;
hold_plots = false;
if length(varargin) > 0
    yrange = varargin{1};
    if ~ isempty(yrange)
        scale_y = true;
    end
end
if length(varargin) > 1
    normalize = varargin{2};
end
if length(varargin) > 2
    hold_plots = varargin{3};
end

%% Get list of file directories:
dpaths = {};
for i = 1:length(evnt);
    pth = evnt(i).dpath;
    tmp = strsplit(pth, '/');
    tmp{end-1} = [tmp{end-1}, '_RAW'];
    out_pth = strjoin(tmp,'/');
    dpaths = [dpaths {out_pth}];
    %% Chang
end
dpaths_set = unique(dpaths);

ecog = [];
for pth = dpaths_set
    block_data = [];
    for i = 1:4
        pth_ch = [pth{1} '/HilbAA_70to150_8band/Wav1' num2str(i) ,'.htk'];
        [ch_data, fs] = readhtk(pth_ch);
        block_data = [block_data; ch_data];
    end
    ecog = [ecog block_data];
end


%% plot based on the relevant flags - works on the assumption of a particular sequence
if hold_plots
    plot_spectrum_chans(ecog, fs, appelation, smoothing, yrange, normalize, hold_plots);
elseif normalize
    plot_spectrum_chans(ecog, fs, appelation, smoothing, yrange, normalize);
elseif scale_y
    plot_spectrum_chans(ecog, fs, appelation, smoothing, yrange);
else
    plot_spectrum_chans(ecog, fs, appelation, smoothing);
end
end
