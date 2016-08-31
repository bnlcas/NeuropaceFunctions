function [Phn] = make_phoneme_erps_rns(evnt, use_HG, sample_rate)
%% This function creates ERPs of the TIMIT sentences for events listed in evnts array
% blocks on the basis of the evnt structure created by DetectEventsQuick

%% Control Phoneme Selection:
Reject_First_Phoneme = true;    % Ignore the ERP from the first phoneme in a sentence (Onset Response)
Reject_Second_Phoneme = true;   % Ignore the Second ERP from the Sentence

%% Initialize ERP Parameters
num_chans = 4;
num_entries = length(evnt);
max_dur_erp = 1.5; % Overestimate of the time duration (s) of TIMIT ERPs
reject_art = true; % Reject events that fall within the interval of a bad time segment.
z_score_prestim = true;

fs = sample_rate;
pre_event_duration = 0.5; % Amount of time to display prior to Onset
ecog_event = zeros(num_chans,max_dur_erp*fs);


Phn.phn_names = [];
Phn.start_times = [];
total_events = num_entries*50; % Overestimate
phn_count = 0;
Phn.ecog = zeros(num_chans,fs*max_dur_erp,total_events);
if reject_art
    reject_erp = false(total_events,1);
end

fs_audio = 16000;

%% Loop Through Events and get start times and directories
start_times = zeros(num_entries,1);          % Event Start Times in Seconds
data_paths = repmat({''},num_entries,1);     % Path of ECoG Data
timit_names = repmat({''},num_entries,1);    % Name of TIMIT File
for i = 1:num_entries
% % %     tmp = strsplit(evnt(i).dpath,'/');
% % %     block_num = tmp{end}
% % %     block_num(1) = [];
% % %     offset = block_offsets(block_num);
    
    start_times(i) = evnt(i).StartTime; %+offset;
    data_paths(i) = {evnt(i).dpath};
    timit_names(i) = {evnt(i).name};
end

%% For Each Phoneme Get a set of Phoneme Names, StartTimes and ECoG Data into Phn
for i = 1:num_entries
    sounddir = '/Users/changlab/Documents/changrepo/matlab/TIMIT/@ECSpeech/Sounds';
    phn_filename = [sounddir, '/',timit_names{i},'.phn'];
    phn_table = import_phn_table(phn_filename);
    
        
    %% Manipulate Onset Phonemes
    if Reject_Second_Phoneme
        phn_table(1:2,:) = [];
    elseif Reject_First_Phoneme
        phn_table(1,:) = [];
    end        
    
    %% Find Start Times
    phn_start_pts = table2array(phn_table(:,1));                        % point in recording where phnome starts
    phn_start_times = ((phn_start_pts/fs_audio) + start_times(i))*1000; % Time in (ms) in audio/ecog where phoneme starts
    Phn.start_times = [Phn.start_times; phn_start_times];
    
    %% Get names
    phn_names = table2array(phn_table(:,3));
    Phn.phn_names = [Phn.phn_names; phn_names];
    
    %% Loop through each phoneme event and Get ECoG ERPs
    for k = 1:length(phn_start_pts)
        phn_count = phn_count + 1;
        if use_HG
            data_path = [data_paths{i}, '/HilbAA_70to150_8band'];
        else
            data_path = data_paths{i};
            data_path_split = strsplit(data_path, '/'); 
            data_path_split{end-1} = [data_path_split{end-1} '_RAW'];
            
            data_path = strjoin(data_path_split,'/');
            data_path = [data_path '/HilbAA_70to150_8band'];
        end
        

        phn_onset_time = phn_start_times(k)-pre_event_duration*1000;
        phn_stop_time = phn_onset_time + max_dur_erp*1000;
        
        if reject_art
            data_path_split = strsplit(data_paths{i}, '/'); 
            data_path_split{end-1} = [data_path_split{end-1} '_RAW'];
            data_path_art = [strjoin(data_path_split,'/') ,'/Artifacts'];
            load([data_path_art, '/badTimeSegments.mat']); % Load the badTimeSegments matrix
            in_interval = @(x,y) (y>= x(1) & y<=x(2));

            for j = 1:size(badTimeSegments,1)
                if (in_interval(badTimeSegments(j,:), phn_onset_time/1000) | in_interval(badTimeSegments(j,:), phn_stop_time/1000))
                    reject_erp(phn_count) = true;
                end
            end
        end

        

        for j = 1:num_chans
            ch_path = [data_path, '/Wav1', num2str(j), '.htk'];
            [ch_data, fs] = readhtk(ch_path,[phn_onset_time, phn_stop_time]);
            dur = length(ch_data);

            ecog_event(j,1:(fs*max_dur_erp)) = ch_data(1:(fs*max_dur_erp)); % an extra point can be added through rounding error
            
        end
        if z_score_prestim
            for j = 1:num_chans
                z_win = [-500 0];
                z_times = (z_win + start_times(i)*1000);
                
                ch_path = [data_path, '/Wav1', num2str(j), '.htk'];
                [ch_data, fs] = readhtk(ch_path,z_times);
                ecog_event(j,:) = (ecog_event(j,:)-mean(ch_data))/std(ch_data);
            end
        end
        
        Phn.ecog(:,:,phn_count) = ecog_event;
    end
end


%% Trim Dead stuff
Phn.ecog(:,:, (phn_count+1):end) = [];
if reject_art
    reject_erp((phn_count+1):end) = [];
    Phn.ecog(:,:,reject_erp) = [];
    Phn.phn_names(reject_erp) = [];
    Phn.start_times(reject_erp) = [];
end

%% Test Plot

test = Phn.ecog(:,:,:);
time_axis = 1000*linspace(-pre_event_duration, max_dur_erp - pre_event_duration, round(max_dur_erp*fs));
Phn.time_axis = time_axis;
% % figure;
% % for i = 1:num_chans
% %     subplot(2,2,i)
% %     shadedErrorBar(time_axis, mean(test(i,:,:),3),nansem(test(i,:,:),3))
% %     xlabel('Time (s)')
% %     ylabel('Awesome Brain Stuff')
% %     title(['Awesome Brain Stuff in Ch ',num2str(i),' For Phoneme ERPs'])
% % end
% % a=1;

