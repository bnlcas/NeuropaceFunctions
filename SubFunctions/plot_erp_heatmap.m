function [] = plot_erp_heatmap(ERPs_ch, varargin)
%% This function draws a heatmap ERP plot for a list of ERPs
% for a single channel,
% a time axis may be provided
% data should be in the format (time_pts x trials) - the most convenient to
% use squeeze in preparation.
% is the ERP structure is provided, ERPs will be sorted
% in Duration from longest to shortest

if isempty(varargin)
    time_axis = 0:(size(ERPs_ch,1)-1);
    order = 1:(size(ERPs_ch,2));
    scale_axis = false;
elseif length(varargin) == 1
    time_axis = varargin{1};
    clip_time = true;
    order = 1:(size(ERPs_ch,2));
    scale_axis = false;
elseif length(varargin) == 2
    time_axis = varargin{1};
    clip_time = true;
    order = 1:(size(ERPs_ch,2));
    ylims = varargin{2};
    if isempty(ylims)
        scale_axis = false;
    else
        scale_axis = true;
    end

else
    time_axis = varargin{1};
    clip_time = true;
    ylims = varargin{2};
    if isempty(ylims)
        scale_axis = false;
    else
        scale_axis = true;
    end
    evnt_struct = varargin{3};
    for i = 1:length(evnt_struct)
        duration(i) = evnt_struct(i).StopTime - evnt_struct(i).StartTime;
    end
    [~,order] = sort(duration);
end
if clip_time
    time_roi = [-500 2000];
    time_inds = [find(min(abs(time_axis-time_roi(1))) == abs(time_axis-time_roi(1)),1), ...
       find(min(abs(time_axis-time_roi(2))) == abs(time_axis-time_roi(2)),1)];
    ERPs_ch = ERPs_ch(time_inds(1):time_inds(2),:);
    time_axis = time_axis(time_inds(1):time_inds(2));
end
    

erps_sort = ERPs_ch(:,order)';
imagesc(flipud(erps_sort))
if scale_axis
    caxis(ylims);
end
ax = gca;
hold on;
[~,zero_ind] = min(abs(time_axis));
plot([zero_ind zero_ind], ax.YLim, 'r')
ylabel('Trial')
xlabel('Time (ms)')

plot_pts = (linspace(1,size(erps_sort,2),6));
%plot_pts = round(plot_pts/10)*10; % round middle points to nearest 10
ax.XTick = plot_pts;
ax.XTickLabel = round(time_axis(round(plot_pts))/5)*5;







