function [x_corr_waves] = audio_env_ecog_corr_grid(evnt, is_hg, zero_wn_pulses,varargin)
%% Plot the X-correlation between the audio envelope and ECoG Data
%% In contrast to the program for the RNS, this program splits each block of ECoG into smaller
% units to take an average amplitude correlation...
% evnt - struct; is the event data structure related to a block of testing
%
% is_hg - boolean, designates wheth
%
% zero_wn_pulses - removes loud pulses of white noise at onset/offset of
% blocks
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

num_blocks = 3;

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
    %% Chop ECoG into smaller blocks
        regions = floor(linspace(1,size(ecog,2), num_blocks+1));
        spacing = min(diff(regions));
        ecog_blocks = [];
        for j = 1:num_blocks
            ecog_blocks = cat(3,ecog_blocks, ecog(:,regions(j):(regions(j)+spacing)));
        end
            
        
    %% Load/Process Audio
        [audio, fs_audio] = readhtk([audio_dir 'ANIN2.htk']);

    %% Chop Audio into smaller blocks
        regions = floor(linspace(1,length(audio), num_blocks+1));       
        spacing = min(diff(regions));
        audio_blocks = [];
        for j = 1:num_blocks
            audio_blocks = [audio_blocks; audio(regions(j):(regions(j)+spacing))];
        end
        

        for j = 1:num_blocks
            % Downsample:
            audio = audio_blocks(j,:);
            ecog = squeeze(ecog_blocks(:,:,j));
            
            fs_ds = fs_ecog;
            audio_ds = resample(audio, 200, 12207);
            if zero_wn_pulses
                pulse_length = 28;
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
                [corrs, lag] = xcorr(ecog(j,:), audio_env, round(fs_ecog*time_window),'coef');
                corr_data = [corr_data; corrs];
            end
            x_corr_waves = cat(3, x_corr_waves, corr_data);
        end
    end

    %% Plot_data
    if plot_output
        tag_title = false;
        scale_y = false;
        if ~isempty(varargin)
            yrange = varargin{1};
            if ~isempty(yrange)
                scale_y = true;
            end
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
