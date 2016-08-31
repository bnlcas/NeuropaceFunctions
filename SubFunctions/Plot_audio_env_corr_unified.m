function [] = Plot_audio_env_corr_unified(varargin)
%% This function takes the ecog audio-envelope cross correlation data
% for several different sources and then plots them on a unified axis
% 
% inputs:
% varagring:
% x_corr_waves a (num_ch x num_timepts x TIMIT_blocks) data matrix
% containing the cross correlation measurement between the ecog and audio
% evelope
% the function accepts a variable number of these, assuming equal number of
% channels and time points
%
% FLAGS:
% 'TimeSpan', the time in seconds on either side of stimulus onset that
% will be plotted
%
% 'YScale', the ylimits.
%
% 'TagTitle', an appelation for the title.
% EX:
% Plot_audio_env_corr_unified(x_corr_waves_1, x_corr_waves_2, 'TimeSpan',1,
% ... 'YScale',[-1.2 1.2], 'TagTitle', 'Broadband');

%% Load Flags
scale_y = false;
timespan = 1;
tag_title = false;
i = 1;
while i<=length(varargin)
    if isstr(varargin{i})
        if strcmpi(varargin{i}, 'YScale')
            scale_y = true;
            y_bounds = varargin{i+1};
            i = i+1; % increment, because next value is known
        end
        if strcmpi(varargin{i}, 'TimeSpan')
            timespan = varargin{i+1};
            i = i+1;
        end
        if strcmpi(varargin{i}, 'TagTitle')
            tag_title = true;
            appelation = varargin{i+1};
            i = i+1;
        end
    end
    i = i + 1;
end

%% Collect and Plot Corr data:
start_flags = false;
i = 1;
color_vals = [0.6 0.2 0.2; 0.2 0.4 0.6; 0.2 0.2 0.8]; % list of colors used for plots
figure; hold on;
while ~start_flags
    if isstr(varargin{i})
        start_flags = true;
    end
    if ~start_flags
        dat = varargin{i};
        means_mat = mean(dat,3);
        sems_mat = nansem(dat,3);
        
        %% Plot data:
        time_axis = 1000*linspace(-timespan,timespan,size(means_mat,2));
        for j = 1:size(means_mat,1)
            subplot(1, size(means_mat,1), j);
            hold on;
            shadedErrorBar(time_axis, means_mat(j,:), sems_mat(j,:), color_vals(:,i),1);
            axis tight;
            plot(get(gca, 'XLim'), [0 0], 'k')
            if scale_y
                ylim(y_bounds)
            end
            plot([0 0], get(gca, 'YLim'), 'k')
            xlabel('Time Offset (ms)')
            ylabel('Correlation')
            if tag_title
                title({[appelation ' Cross-Correlation between']; ['Audio Evelope & ECoG in Ch ', num2str(j)]})
            else
                title({['Cross-Correlation between Audio']; ['Evelope & ECoG in Ch ', num2str(j)]})
            end
        end
  
    end

        
        
        
    i = i+1;
    if i>length(varargin)
        start_flags = true;
    end
end






end
