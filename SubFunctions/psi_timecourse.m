function [out_matrix] = psi_timecourse(phn, appelation)
%% take the phoneme event structure and plots the time courses for the PSI by sliding the window

%% parameters:
fs = 250; % 250 hz sample rate
fs = round(1000*length(phn.time_axis)/(max(phn.time_axis)-min(phn.time_axis)));

num_phonemes = 36; % Hard coded magic constant - reference to list found in PSI_RNS
window_width = 0.05;
tmin = 0;
tmax = 0.5;
times_center = linspace(tmin, tmax, round(fs*(tmax-tmin)/5));

%% Get PSI's
out_matrix = zeros(length(times_center), num_phonemes, size(phn.ecog,1));
for i = 1:length(times_center)
    [tmp, phn_names] = PSI_RNS(phn, fs, false, times_center(i), window_width);
    out_matrix(i,:,:) = tmp;
end

%% Plot data:
plot_chans = 4;
figure;
for i = 1:plot_chans
    subplot(1,4,i)
    dat = squeeze(out_matrix(:,:,i))'; % time_pts x phoneme
    imagesc(dat)
    ax = gca;
    xpts = round(linspace(1,length(times_center),5));
    ax.XTick = xpts;
    ax.XTickLabels = 1000*times_center(xpts);
    xlabel('time (ms)')
    
    ax.YTick = 1:length(phn_names);
    ax.YTickLabel = phn_names;
    title([appelation ' PSI timecourse for Ch ' num2str(i)]);
end
a=1;

        