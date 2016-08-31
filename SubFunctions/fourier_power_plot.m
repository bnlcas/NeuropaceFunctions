function [] = fourier_power_plot(signal, sample_rate, varargin);
%% Generates a plot of the power distribution of a input signal
%Plot the magetude of the frequency components up to the nyquist frequency
%
% INPUTS:
%
% Signal - data to be transformed
% 
% sample_rate
%
% Variable Inputs:
%
% 1 - smoothing factor (amount of smoothing of spectrogram
%
% 2 - normalize to frequency band - setting to normalize data to the power for a given
% frequecy band, defaults to no normalization if empty

smoothing_param = 500; % Default setting
normalize = false;
if length(varargin) > 0
    smoothing_param = varargin{1};
end
if length(varargin) > 1
    normalization_window = varargin{2};
    if ~isempty(normalization_window)
        normalize = true;
    end
end

[ps_density,freq_axis] = pwelch(signal,smoothing_param,[],[], sample_rate); %smoothing_param/2,smoothing_param,sample_rate);

%freq = fft(signal);
%sample_freqs = floor(length(signal)/2);
%nyquist = sample_rate/2;
%freq_axis = linspace(0,nyquist,sample_freqs);

%ps_density = freq.*conj(freq);

if normalize
    [~,roi(1)] = min(abs(freq_axis - normalization_window(1)));
    [~,roi(2)] = min(abs(freq_axis - normalization_window(2)));
    power_per_bin = mean(ps_density(roi(1):roi(2)));
    ps_density = ps_density/power_per_bin;
end

% smooth
%ps_density_smth = smooth(ps_density, smoothing_param);

%semilogy(freq_axis, ps_density_smth(1:sample_freqs));
semilogy(freq_axis, ps_density);
loglog(freq_axis, ps_density);
xlim([0 sample_rate/2])
%plot(freq_axis, 1./ps_density_smth(1:sample_freqs));

xlabel('Frequency (Hz)')
ylabel('Power')
if normalize
    ylabel('Normalized Power')
end
title('Spectrum of Incident Signal')

%% Additional code to confirm ~ 1/f dependence
% % %ylabel('Inverse Power Density (1/P)')
% % hold on;
% % freq_bounds = [50 100];
% % in_bounds  = freq_axis>freq_bounds(1) & freq_axis<=freq_bounds(2);
% % x_dat = [ones(sum(in_bounds),1), freq_axis(in_bounds)'];
% % y_dat = 1./ps_density_smth(1:sample_freqs);
% % y_dat = y_dat(in_bounds);
% % coeffs = regress(y_dat, x_dat);
% % plot(freq_bounds, [(freq_bounds(1)*coeffs(2)+coeffs(1)), (freq_bounds(2)*coeffs(2)+coeffs(1))],'r','LineWidth',2)
% % ylim([0 1.5*(freq_bounds(2)*coeffs(2)+coeffs(1))])

end
