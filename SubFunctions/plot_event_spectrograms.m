function [] = plot_event_spectrograms(erps_bb, time_axis);
%This function takes the broadband ecog data from a ch x time x trials data matrix
% and generates a spectrogram for each channel over time by averaging
% across trials

fs = length(time_axis)/(max(time_axis)/1000 - min(time_axis)/1000);

figure;
for i = 1:size(erps_bb,1)
    ch_dat = squeeze(erps_bb(i,:,:));
    event_spectrograms = [];
    for j = 1:size(ch_dat,2)
        [s, f, t] = spectrogram(ch_dat(:,j), 128,120,128,fs,'yaxis');
        event_spectrograms = cat(3, event_spectrograms, s);
    end
    mean_spectrogram = mean(abs(event_spectrograms),3);
    % log scale data (used in matlab's defaul)
    log_spectrum = log(mean_spectrogram);
    %log_spectrum = gdivide(abs(mean_spectrogram), abs(sum(mean_spectrogram,2)));
    subplot(1,size(erps_bb,1),i)
    
    log_spectrum = flipud(log_spectrum);
    f = sort(f, 'descend');
    t_span = linspace(min(time_axis), max(time_axis), length(t));
    imagesc(t_span, f, log_spectrum)
    set(gca,'YDir','normal')

    ylim([0 125]);
    
    %% Plot time of onset
    hold on;
    plot([0 0], get(gca, 'YLim'), 'r')
    xlim([-500 1000])
    xlabel('Time (ms)')
    ylabel('Frequency (hz)')
    title(['Spectral ERPs in Ch ' num2str(i)])
end
