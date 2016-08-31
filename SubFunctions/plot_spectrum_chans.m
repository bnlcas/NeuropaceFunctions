function [] = plot_spectrum_chans(raw_data, sample_rate, title_description, varargin)
%% Plots a 1 x 4 set of Spectrograms for raw data:
%
% INPUTS:
%
% raw_data - raw data for spectrugrams
% sample_rate
% title_description - added to plot titles
% Varargin order:
%

% 1: smoothing level - amount of spectral smoothing
% 2: scale y axis - specify the scale of the yaxis
% 3: normalize - boolean to normalize spectrum to power
% 4: hold_plots - boolean to overwrite existing figure or plot on a new
% figure

if isempty(title_description)
    title_description = 'ENTER DESCRIPTION';
end

num_plots = 4;
scale_y = false;
normalize = false;
hold_plots = false;
if isempty(varargin)
    smoothing_level = 100;
elseif length(varargin) > 0
    smoothing_level = varargin{1};
end
    
if length(varargin) > 1
    yrange = varargin{2};
    if ~isempty(yrange)
        scale_y = true;
    end
end
if length(varargin) > 2
    normalize = varargin{3};
end
if length(varargin) > 3
    hold_plots = varargin{4};
end

if ~hold_plots
    figure;
end

for j = 1:num_plots
     subplot(1, num_plots, j)
     if hold_plots
         hold on;
     end
     if normalize
        normalize_freqs = [0, 125];
     else
         normalize_freqs = [];
     end
     %ecog_ch = squeeze(raw_data(j,:,:));
     %ecog_ch = ecog_ch(:);

     ecog_ch = raw_data(j,:);
   
     fourier_power_plot(ecog_ch, sample_rate, smoothing_level, normalize_freqs);
     set(gca, 'XLim', [0 125])
     if scale_y
         set(gca,'YLim', yrange)
     end
     title(['Spectrum of ' title_description ' in CH ' num2str(j)])
end