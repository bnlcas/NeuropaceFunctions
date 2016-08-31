function [ecog_speech, ecog_nonspeech] = gather_raw_speech_nonspeech_ecog(evnt)
% This function takes an evnt strucutre for a block of testing and parses
% the raw data to create two num_chans x timepts data matrices (one for speech events
% and one for non-speech times) that can be used for later processing.



num_chans = 4;

%% Get list of block directories:
block_dirs = cell(length(evnt),1);
for i = 1:length(evnt)
    block_dirs{i} = evnt(i).dpath;
end
block_dirs = unique(block_dirs);


%% For each block assemble a matrix of non_speech and speech data;
ecog_speech = [];
ecog_nonspeech = [];
for i = 1:length(block_dirs)
    start_times = [];
    stop_times = [];
    for j = 1:length(evnt)
        % If in block
        if strcmp(evnt(j).dpath, block_dirs(i))
            % Reformat dpath to RAW dir
            dpath = strsplit(evnt(j).dpath,'/');
            dpath{end-1} = [dpath{end-1} '_RAW'];
            dat_dir = [strjoin(dpath,'/'), '/HilbAA_70to150_8band/'];
            
            
            start_times = [start_times; evnt(j).StartTime];
            stop_times = [stop_times; evnt(j).StopTime];
        end
    end
    start_times_speech = sort(start_times);
    stop_times_speech = sort(stop_times);
    
    %% Load data and compliment
    ecog_block = [];
    for j = 1:num_chans
        [dat,fs] = readhtk([dat_dir ,'Wav1' num2str(j) '.htk']); % Technically relies on memory effects...
        ecog_block = [ecog_block;dat];
    end
        
    start_times_nonspeech = [0; stop_times_speech];                 % Assumes Block starts with no speech
    stop_times_nonspeech = [start_times_speech; (length(dat)/fs)];  % Assumes Block ends with no speech

    %% Block data:
    for j=1:length(start_times_speech)
        ecog_speech = [ecog_speech, ecog_block(:,round(fs*start_times_speech(j)):round(fs*stop_times_speech(j)))];
    end
    for j = 1:length(start_times_nonspeech)
        ecog_nonspeech = [ecog_nonspeech, ecog_block(:,(1+round(fs*start_times_nonspeech(j))):round(fs*stop_times_nonspeech(j)))];
    end
end




end

