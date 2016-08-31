function erps_zscore = z_score_erps_prestim(erps, time_axis)
% z-normalizes erps to the time window of -0.5 to 0 seconds prior to the
% stimulus onset

%% Get ROI:
z_roi = [-0.5 0];
[~,z_inds(1)] = min(abs(time_axis-z_roi(1)));
[~,z_inds(2)] = min(abs(time_axis-z_roi(2)));

erps_zscore = zeros(size(erps));
% for each channel and each trial, z-score data:
for i = 1:size(erps,1)
    for j = 1:size(erps,3)
        data = squeeze(erps(i,:,j));
        z_data = data(z_inds(1):z_inds(2));
        erps_zscore(i,:,j) = (data - mean(z_data))/std(z_data);
    end
end


end
