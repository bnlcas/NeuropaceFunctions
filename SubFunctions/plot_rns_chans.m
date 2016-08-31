function [] = plot_rns_chans(erps, time_axis, is_hg, title_description, yrange)
%% Plots a 1 x 4 set of ERPs for RNS data:
if isempty(title_description)
    title_description = 'ENTER DESCRIPTION';
end
    
num_plots = 4;
time_range = [-200 1000];
clip_data = true; % parameter to clip data to scaling window
if clip_data
    twin(1) = find(time_axis <= time_range(1),1,'last');
    twin(2) = find(time_axis >= time_range(2),1);
    erps = erps(:,twin(1):twin(2),:);
    time_axis = time_axis(twin(1):twin(2));
end
    


%% Plot channels
figure;

for j = 1:num_plots
     subplot(1, num_plots, j)
      ecog_ch = squeeze(erps(j,:,:));
     
     shadedErrorBar(time_axis, mean(ecog_ch,2), nansem(ecog_ch,2));
     axis tight
    
%     if exist('yrange', 'var')
     if ~isempty(yrange)
        ylim(yrange);
     end
     
     if is_hg
         title([title_description ' High Gamma in CH ' num2str(j)])
         ylabel('Z-normalized High Gamma Intensity')
     else
         title([title_description ' Broadband in CH ' num2str(j)])
         ylabel('Z-normalized Broadband Intensity')
     end
     xlabel('Time (ms)')
     hold on;
     ax = gca;
     plot(ax.XLim, [0 0],'k')
     plot([0 0], ax.YLim, 'k')
     xlim(time_range)
 end