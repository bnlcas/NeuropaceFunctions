function [x_corr_waves] = audio_env_ecog_corr(evnt, is_hg, zero_wn_pulses)
%% Plot the X-correlation between the audio envelope and ECoG Data
%
% evnt - struct; is the event data structure related to a block of testing
%
% is_hg - boolean, designates wheth
%
% zero_wn_pulses - removes loud pulses of white noise at onset/offset of
% blocks
%% System Paramters:
time_window = 1; % 1 second of cross correlation.
num_chans = 4;
plot_output = true;

%% Load directory data
dpaths = {};
for i = 1:length(evnt)
    dpaths{i} = evnt(i).dpath;
end
dpaths = unique(dpaths);


%% Load Data and calcualte Xcorrs for each block
x_corr_waves = [];
for i = 1:length(dpaths)
    pth = dpaths{i};
    audio_dir = [pth, '/Analog/']; % Analog is only in the HG (Non-RAW dir)
    if ~is_hg
        tmp = strsplit(pth, '/');
        tmp{end-1} = [tmp{end-1}, '_RAW'];
        pth = strjoin(tmp,'/');
    end
    ecog_dir = [pth, '/HilbAA_70to150_8band/'];
    
    %% Load ECoG
        ecog = [];
        for j = 1:num_chans
            [data, fs_ecog] = readhtk([ecog_dir 'Wav1' num2str(j) '.htk']);
            ecog = [ecog; data];
        end
        
        
    %% Load/Process Audio
        [audio, fs_audio] = readhtk('ANIN2.htk');

        % Downsample:
        fs_ds = fs_ecog;
        audio_ds = resample(audio, fs_ds, fs_audio);
        if zero_wn_pulses
            pulse_length = 27;
            pulses_on = false(1,length(audio_ds));
            pulses_on(1:(fs_ds*pulse_length)) = true;
            pulses_on((end-(fs_ds*pulse_length)):end) = true;
        else
            pulses_on = false(1,length(audio_ds));
        end
        audio_ds(pulses_on) = 0;
        audio_env = abs(hilbert(audio_ds));
        
        ecog(:,pulses_on) = 0;
     %% Trim data to make length_equivalent
        data_size = min(length(audio_env),size(ecog,2));
        ecog = ecog(:,1:data_size);
        audio_env = audio_env(1:data_size);
        pulses_on = pulses_on(1:data_size);
        
        %Clip Pulse data:
        ecog(:,pulses_on) = [];
        audio_env(pulses_on) = [];
     
     %% Calc x_corr:
        corr_data = [];
        for j = 1:num_chans
            [corrs, lag] = xcorr(ecog(j,:), audio_env, fs_ecog*time_window,'coef');
            corr_data = [corr_data; corrs];
        end
        x_corr_waves = cat(x_corr_waves, corr_data, 3);
end

%% Plot_data
if plot_output
    time_axis = linspace(-fs_ecog*time_window,fs_ecog*time_window,size(x_corr_waves,2));
    figure;
    for j = 1:num_chans
        subplot(1, num_chans, j)
        shadedErrorBar(time_axis, mean(x_corr_waves,3), nansem(x_corr_waves,3))
    end
end



end
