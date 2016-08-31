function [] = multi_plot_rns_chans(title_description, yrange, varargin)
%% Plots a 1 x 4 set of ERPs for RNS data:
% inputs
%
% title_descrition - title for each subplot
%
% varargin:
% enter alternating (erps, time_axis), inputs
normalize = false;

if mod(length(varargin),2) ~= 0
    throw('Enter pairings of time axis and ecog data jerkface')
else
    num_data_sources = length(varargin)/2;
end
if isempty(yrange)
    yrange = [-1 1];
end
if isempty(title_description)
    title_description = 'ENTER DESCRIPTION';
end

for j = 1:num_data_sources
    ecog_ind = 2*j-1;
    time_axis_ind = 2*j;
    
    ecog{j} = varargin{ecog_ind};
    time_axis{j} = varargin{time_axis_ind};
    
end

if normalize
    for i = 1:num_data_sources
        ecog_means = mean(ecog{i},3);
        max_mean_ecog = max(abs(ecog_means(:)));
        %max_mean_ecog = max(abs(ecog_means,[]2);
        %max_mean_ecog = max(abs(ecog_means),[],2);
        ecog{i} = ecog{i}/max_mean_ecog;
        %ecog{i} = gdivide(ecog{i},max_mean_ecog);
    end
end
    
num_plots = 4;
time_range = [-200 1000];


%% Clip data:
clip_data = true; % parameter to clip data to scaling window
if clip_data
    for j = 1:num_data_sources    
        ecog_tmp = ecog{j};
        time_axis_tmp = time_axis{j};
        twin(1) = find(time_axis_tmp <= time_range(1),1,'last');
        twin(2) = find(time_axis_tmp >= time_range(2),1);
        
        ecog{j} = ecog_tmp(:, twin(1):twin(2) ,:);
        time_axis{j} = time_axis_tmp(twin(1):twin(2));
    end
end



%% Plot Data
figure;
if num_plots == 1
    colorlist = [0 0 0];
else
    colorlist = [0.4 0.4 0.2; 0.8 0.1 0; 0 0.1 0.8; 0.1 0.1 0.1];
end
for j = 1:num_plots
     subplot(1, num_plots, j)
     hold on;
     for i = 1:num_data_sources
         tmp = ecog{i};
         ecog_ch = squeeze(tmp(j,:,:));
         time_axis_plot = time_axis{i};
     
         shadedErrorBar(time_axis_plot, mean(ecog_ch,2), nansem(ecog_ch,2), colorlist(i,:), 1);
         axis tight
   
         if ~isempty(yrange)
            ylim(yrange);
         end

         title({title_description; ['in Ch ' num2str(j)]})
         ylabel('Normalized Intensity')
         xlabel('Time (ms)')
         ax = gca;
         plot(ax.XLim, [0 0],'k')
        plot([0 0], ax.YLim, 'k')
        xlim(time_range)
     end
 end
