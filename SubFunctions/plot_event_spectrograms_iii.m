function [] = plot_event_spectrograms_iii(erps_bb, time_axis);
% This function takes the broadband ecog erp data from a (Channel x Time x Trials) data matrix
% along with the time axis of that data and generates a spectrogram for
% each of the channels (default n = 4).
% Since the hilbert transform is known to have edge effects, windowing, zscoring, or
% clipping data is recommended. In practive, using a longer time axis erp
% was selected since this allows the edge effects to be safely ignored.
%
% The spectrogram is formed by the following method: for logarithmically spaced frequencies,
% the raw data is band pass filterd, and the analytic amplitude of this
% band is then calculated - this is used to plot the intensity of the ECoG
% at this frequency over time. Note that the analytic amplitude of each
% such band is z-scored according to the band's activity in the 500 ms
% prior to the onset of the stimulus.
%
% inputs:
% erps_bb - channels x time_points x trials matrix of broadband ecog data
%
% time_axis - 1xtime_points data listing the time (relative to event onset)
% of each time point in the erps_bb input
%
% variable parameters (should in incorporated in function call once a stable
% release has been reached)
%
% color axis scaling (line 106)
% current setting: caxis([0 1]) - (feature recently added during debugging)


%% Set frequency band parameters (copy and pasted in).
fs_ecog = length(time_axis)/(max(time_axis)/1000 - min(time_axis)/1000);
freq_axis = (fs_ecog/2)*(0:floor(size(erps_bb,2)/2-1))/floor(size(erps_bb,2)/2-1);
freq_axis = [freq_axis, sort(freq_axis,'descend')];


freq_samples = [0.0180000000000000,0.0703239906734951,0.173746807215294,0.336310300587718,0.562480141016083,0.854974977634691,1.21558791985701,1.64557736462310,2.14586959638287,2.71717228814078,3.36004234112678,4.07492865382648,4.49908599565877,4.96739366892354,5.48444726014817,6.05532070822666,6.68561609588494,7.38151862391541,8.14985730765295,8.99817199131751,9.93478733784706,10.9688945202963,12.1106414164533,13.3712321917698,14.7630372478308,16.2997146153059,17.9963439826350,19.8695746756941,21.9377890405926,24.2212828329065,26.7424643835396,29.5260744956615,32.5994292306116,35.9926879652699,39.7391493513881,43.8755780811850,48.4425656658129,53.4849287670791,59.0521489913229,65.1988584612231,71.9853759305396,79.4782987027759,87.7511561623698,96.8851313316256,106.969857534158,118.104297982645];
window_widths = [0.0523239906734951,0.103422816541799,0.162563493372424,0.226169840428365,0.292494836618608,0.360612942222318,0.429989444766091,0.500292231759772,0.571302691757911,0.642870052985993,0.714886312699707,0.787271648319059,0.827230910894714,0.869218371321770,0.913336974105689,0.959694888868996,1.00840577556066,1.05958906312661,1.11337024232463,1.16988117340155,1.22926040938710,1.29165353579707,1.35721352757867,1.42610112417325,1.49848522361586,1.57454329663811,1.65446182178942,1.73843674264354,1.82667394821137,1.91938977773799,2.01681155112132,2.11917812625322,2.22674048464926,2.33976234680310,2.45852081877419,2.58330707159413,2.71442705515734,2.85220224834648,2.99697044723171,3.14908659327622,3.30892364357884,3.47687348528707,3.65334789642274,3.83877955547597,4.03362310224263,4.23835625250643];

freq_samples(1:7) = []; % data should contain 39 bands(standard 40 with one dropped because fs = 250)
window_widths(1:7) = [];

% make data an even number of points
if mod(size(erps_bb,2),2) == 1
    erps_bb(:,end,:) = [];
end


% MAKE BANDS
bands_data = zeros(size(erps_bb,1),size(erps_bb,2), size(erps_bb,3),length(freq_samples));
for k = 1:length(freq_samples)
    %% Filter ECoG Channels
    gaussian_window = @(f) exp(-((f-freq_samples(k)).^2)/(2*window_widths(k).^2));
    for i = 1:size(erps_bb,1)
        for j = 1:size(erps_bb,3)
            %ecog_ch = erp_win.*squeeze(erps_bb(i,:,j));
            ecog_ch = squeeze(erps_bb(i,:,j));
            
            % normalization bounds (-500 to 0 ms prestimulus).
            [~,t1] = min(abs(time_axis- (-500))); %
            [~,t2] = min(abs(time_axis));
            ecog_ch = (ecog_ch-mean(ecog_ch(t1:t2)))/std(ecog_ch(t1:t2));
            
            ecog_fft = fft(ecog_ch);

            T=length(ecog_ch);
            h = zeros(1,T);
            if 2*fix(T/2)==T %if T is even
                h([1 T/2+1]) = 1;
                h(2:T/2) = 2;
            else
                h(1) = 1; h(2:(T+1)/2) = 2;
            end

            gauss_filter = gaussian_window(freq_axis);

            hilbdata=ifft(ecog_fft(end,:).*(gauss_filter.*h),T);
            envData=abs(hilbdata);

            bands_data(i,:,j,k) = (envData - mean(envData(t1:t2)))/std(envData(t1:t2)); % Z-score band prestimulus
            %bands_data(i,:,j,k) = envData; 
        end
    end


end

%% AVERGAGE AND PLOT
spectral_erps = squeeze(mean(bands_data,3));

%% Clip Time edges (necessary because of edge effects of hilbert transform
edge_clip = floor(length(time_axis)*0.1/2);
spectral_erps(:,1:edge_clip,:) = [];
spectral_erps(:,(end-edge_clip):end,:) =[];
time_axis(1:edge_clip) = [];
time_axis((end-edge_clip):end) = [];

[~,t1] = min(abs(time_axis- (-200))); % time bounds for display (necessary for use of imagesc
[~,t2] = min(abs(time_axis -1000));

figure;
for i = 1:size(spectral_erps,1)
    subplot(1,4,i)
    ch_spectrogram = squeeze(spectral_erps(i,:,:));

    ch_spectrogram = flipud(ch_spectrogram'); % flip data for xy axis to be time/freq
    imagesc(ch_spectrogram)
    caxis([0 1]) % color limitation (used in debug phases)
    ax = gca;

    xlim([t1, t2])
    
    % set x ticks and y tick labels
    xticks = round(linspace(t1,t2,5));
    for k = 1:length(xticks)
        xlabs{k} = num2str(floor(time_axis(xticks(k))));
    end
    ax.XTick = round(xticks);
    ax.XTickLabel = xlabs;
    
    
    yticks = 1:6:length(freq_samples);
    ax.YTick = yticks;
    for k = 1:length(yticks)
        ylabs{k} = num2str(round(freq_samples(1+length(freq_samples)-yticks(k))));
    end
    ax.YTickLabel = ylabs;
    
    % plot line marking event onset time
   hold on;
   [~,z_ind] = min(abs(time_axis));
   plot([z_ind, z_ind], ax.YLim,'r')
    
   xlabel('time (s)')
    ylabel('Frequency (hz)')
    title(['Z Scored Intensity of Spectral Bands of Ch ' num2str(i)])
end




end

