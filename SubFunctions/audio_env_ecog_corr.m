function [x_corr_waves] = audio_env_ecog_corr(evnt, is_hg, zero_wn_pulses,varargin)
%% Plot the X-correlation between the audio envelope and ECoG Data
%
% evnt - struct; is the event data structure related to a block of testing
%
% is_hg - boolean, designates wheth
%
% zero_wn_pulses - removes loud pulses of white noise at onset/offset of
% blocks (acutually deletes them)
%
% varargin:
% 
% 1 - ylim
% 
% 2 - plot title label
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
            if ~is_hg
                data = abs(hilbert(data)); % Take the analtitic envelope of broad band data
                data = (data - mean(data))/std(data);
            end
            ecog = [ecog; data];
        end
        
        
    %% Load/Process Audio
        [audio, fs_audio] = readhtk([audio_dir 'ANIN2.htk']);

    %% Clear Whitenoise Pluses
        if zero_wn_pulses
            pulse_length = 28;
            pulses_on_audio = false(1,length(audio));
            pulses_on_audio(1:(fs_audio*pulse_length)) = true;
            pulses_on_audio((end-(fs_audio*pulse_length)):end) = true;
            
            pulses_on_ecog = false(1,size(ecog,2));
            pulses_on_ecog(1:(fs_ecog*pulse_length)) = true;
            pulses_on_ecog((end-(fs_ecog*pulse_length)):end) = true;
        else
            pulses_on_audio = false(1,length(audio_ds));
            pulses_on_ecog = false(1,size(ecog,2));
        end
        audio(pulses_on_audio) = [];
        ecog(:,pulses_on_ecog) = [];


    %% Get Audio Envelope
    % Downsample:
        fs_ds = fs_ecog;
        %audio_env_ds1 = abs(hilbert(resample(audio, fs_ds, fs_audio)));
        %audio_env_ds2 = resample(abs(hilbert(audio)), fs_ds, fs_audio);

        %audio_env = audio_env_ds1;
        %audio_env = audio_env_ds2;
        audio_env = resample(abs(hilbert(audio)), fs_ds, fs_audio);
        
        
        %% Scientific Seqway
% % % % %         audio_ds = resample(audio, fs_ds, fs_audio);
% % % % %         time_axis_ds = linspace(0,length(audio_ds)/fs_ds, length(audio_ds));
% % % % % 
% % % % %         figure; subplot(1,3,1)
% % % % %         plot(time_axis, audio, time_axis_ds, audio_ds); xlabel('time (s)'), ylabel('Intensity')
% % % % %         legend('audio', 'downsampled')
% % % % %         xlim([51 51.5])
% % % % %         title('Comparison of Audio and Downsampled Audio')
% % % % %         
% % % % %         subplot(1,3,2)
% % % % %         plot(time_axis, audio, time_axis_ds, audio_env_ds2) %, time_axis_ds, audio_env_ds1);
% % % % %          xlabel('time (s)'); ylabel('Intensity')
% % % % %          xlim([51 51.5])
% % % % %          title('Comparison of Audio and Downsampled Envelope')
% % % % %          legend('Audio', 'Envelope')
% % % % %          
% % % % %               subplot(1,3,3)
% % % % %         plot(time_axis, audio, time_axis_ds, audio_env_ds1) %, time_axis_ds, audio_env_ds1);
% % % % %          xlabel('time (s)'); ylabel('Intensity')
% % % % %          xlim([51 51.5])
% % % % %          title('Comparison of Audio and Evnvelope of Downsampled Audio')
% % % % %          legend('Audio', 'Envelope')
% % % % %          


        
        
     %% Trim data to make length_equivalent
        data_size = min(length(audio_env),size(ecog,2));
        ecog = ecog(:,1:data_size);
        audio_env = audio_env(1:data_size);
 %       pulses_on = pulses_on(1:data_size);
        
% % % %         %Clip Pulse data:
% % % %         ecog(:,pulses_on) = [];
% % % %         audio_env(pulses_on) = [];
     
     %% Calc x_corr:
        corr_data = [];
        for j = 1:num_chans
            [corrs, lag] = xcorr(ecog(j,:), audio_env, round(fs_ecog*time_window),'coef');
            corr_data = [corr_data; corrs];
        end
        x_corr_waves = cat(3, x_corr_waves, corr_data);
end

%% Plot_data
if plot_output
    tag_title = false;
    scale_y = false;
    if ~isempty(varargin)
        scale_y = true;
        yrange = varargin{1};
        if length(varargin) > 1
            tag_title = true;
            appelation = varargin{2};
        end
    end
    time_axis = 1000*linspace(-time_window,time_window,size(x_corr_waves,2));
    figure;
    for j = 1:num_chans
        subplot(1, num_chans, j)
        shadedErrorBar(time_axis, mean(squeeze(x_corr_waves(j,:,:)),2), nansem(squeeze(x_corr_waves(j,:,:)),2));
        axis tight;
        hold on;
        plot(get(gca, 'XLim'), [0 0], 'k')
        if scale_y
            ylim(yrange)
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



end
