function [ERPs, time_axis] = make_sentence_erps_vns(evnt, use_HG, varargin)
%This function creates ERPs of the TIMIT sentences for each of the ten
%blocks on the basis of the evnt structure created by DetectEventsQuick
%
% Inputs:
% use_HG - determines whether to use the Raw/High Gamma data (both must be
% formatted to RNS standard
%
%
% variable input:
%
% 1.)fs - sample frequency of data (relevant 
% 2.) time_winow (in seconds)

%% function Parameters:
TIMIT_block_delay = 0.12;
num_chans = 4;
num_entries = length(evnt);
max_length_ecog = 2; % Overestimate of the duration of TIMIT ERPs
max_dur = 2; % Test for how long data actually is
fs = 250;

pre_event_duration = 1; % Amount of time to display prior to Onset
post_event_delay = 1; % Amount of Time to display after end
reject_art = true; % Reject events that fall within the interval of a bad time segment.
zscore_pretrial = true; % zscore data to the -0.5 to 1 interval prior to onset

if length(varargin) > 0
    fs = varargin{1};
end
if length(varargin) > 1
    twin = varargin{2};
    pre_event_duration = -twin(1);
    post_event_delay = twin(2);
end

ERPs = zeros(num_chans, fs*(post_event_delay + pre_event_duration), num_entries);
if reject_art
    reject_erp = false(num_entries,1);
end


for i = 1:num_entries
    evnt(i).StartTime = evnt(i).StartTime + TIMIT_block_delay;
    start_time = (evnt(i).StartTime - pre_event_duration)*1000;
    start_times(i) = evnt(i).StartTime;
    %stop_time = (evnt(i).StopTime + post_event_delay)*1000;
    stop_time = (evnt(i).StartTime+post_event_delay)*1000;
    data_path = evnt(i).dpath;
    %% Determine whether artifacts are rejected
    if reject_art
        data_path_split = strsplit(data_path, '/'); 
        data_path_split{end-1} = [data_path_split{end-1} '_RAW'];
        data_path_art = [strjoin(data_path_split,'/') ,'/Artifacts'];
        load([data_path_art, '/badTimeSegments.mat']); % Load the badTimeSegments matrix
        in_interval = @(x,y) (y>= x(1) & y<=x(2));
        
        for j = 1:size(badTimeSegments,1)
            if (in_interval(badTimeSegments(j,:), start_time/1000) | in_interval(badTimeSegments(j,:), stop_time/1000))
                reject_erp(i) = true;
            end
        end
    end

            
    if use_HG
        data_path = [data_path, '/HilbAA_70to150_8band'];
    else
        data_path_split = strsplit(data_path, '/'); 
        data_path_split{end-1} = [data_path_split{end-1} '_RAW'];

        data_path = strjoin(data_path_split,'/');
        data_path = [data_path '/HilbAA_70to150_8band'];
    end
    
    for j = 1:num_chans
        ch_path = [data_path, '/Wav1', num2str(j), '.htk'];
        [ch_data, fs] = readhtk(ch_path,[start_time, stop_time]);
        dur = length(ch_data);
        

        ERPs(j,1:dur,i) = ch_data;
        
%         if dur > max_dur
%             dur = max_dur;
%         end
    end
    
end

%% Trim Trailing Zeros
ERPs(:,(dur+1):end,:) = [];
test = mean(ERPs(:,:,51:100),3);
time_axis = 1000*(linspace(-pre_event_duration, post_event_delay, dur));
%figure; plot(repmat(time_axis,4,1)',test')
a=1;


%% Eliminate Bad Trials:
if reject_art
    ERPs(:,:,reject_erp) = [];
end

if zscore_pretrial
    %erps_zscore = z_score_erps_prestim(erps, time_axis)
    z_bounds = [-500 0];
    z_time_range = (time_axis >= z_bounds(1)) & (time_axis <= z_bounds(2));
    for i = 1:size(ERPs,1)
        for j = 1:size(ERPs,3)
            dat = squeeze(ERPs(i,:,j));
            %z_dat = ERPs(:,z_time_range,:);
            ERPs(i,:,j) = (dat - mean(dat(z_time_range)))/std(dat(z_time_range));
        end
    end
end