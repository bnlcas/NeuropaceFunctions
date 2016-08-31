function [r, lag] = xcorr_all_data(use_HG)
%% Takes a n chann x time_pts x trials object and zsocres
basepath = '/Users/changlab/Documents/data/EC71';
if use_HG
    fs = 400;
else
    fs = 254;
end

channels = [41 33 25 17]; % Channels of interest
ref_chans = [49 41 33 25];
%ref_chans = [42 34 26 18];
ref_chans = [];
chans = [channels, ref_chans];

%% Load Audio:
audio_path = [basepath, '/EC71_B1/Analog/ANIN1.htk'];   % Audio 1 is a More accurate reflection of what patient Heard
[audio_b1, fs_audio] = readhtk(audio_path);
audio_path = [basepath, '/EC71_B5/Analog/ANIN1.htk'];
[audio_b2, fs_audio] = readhtk(audio_path);
audio_data =  [audio_b1, audio_b2];
audio_data = smooth(abs(hilbert(audio_data)),fs_audio/100); %10 ms window
if use_HG
    audio_data = resample(audio_data,fs, round(fs_audio));
else
    audio_data = decimate(audio_data, round(fs_audio/fs));
end
audio_data = (audio_data-mean(audio_data))/std(audio_data);


%% Load ECoG
ecog = [];
for i = 1:length(chans)
    if use_HG
        data_path_B1 = '/Users/changlab/Documents/data/EC71/EC71_B1/HilbAA_70to150_8band';
        data_path_B5 = '/Users/changlab/Documents/data/EC71/EC71_B5/HilbAA_70to150_8band';
    else
        data_path_B1 = '/Users/changlab/Documents/data/EC71/EC71_B1/RawHTK';
        data_path_B5 = '/Users/changlab/Documents/data/EC71/EC71_B5/RawHTK';
    end
    
    if use_HG
        data1 = mean(readhtk([data_path_B1, '/Wav2', num2str(chans(i)),'.htk']),1);
        data5 = mean(readhtk([data_path_B5, '/Wav2', num2str(chans(i)),'.htk']),1);
        data = [data1, data5];
    else
        [data1, fs_ecog] = readhtk([data_path_B1, '/Wav2', num2str(chans(i)),'.htk']);
        data5 = readhtk([data_path_B5, '/Wav2', num2str(chans(i)),'.htk']);    
        data = [data1, data5];
        data = smooth(abs(hilbert(data)),fs_ecog/100); %10 ms window
        data = decimate(data,12);
    end
    %Zscore:
    data = (data-mean(data))/std(data);
    ecog_mp(i,:) = data;
end
%% Bipolar Junction on ECoG
%ecog = ecog_mp(1:(length(chans)/2),:) - ecog_mp((1+length(chans)/2):end,:);
ecog = ecog_mp

%% Clip DownSampling Errors
DS_Thresh = 2;
if(length(audio_data) ~= size(ecog,2))
    if (abs(length(audio_data) - size(ecog,2)) <=DS_Thresh)
        max_size = min(length(audio_data), size(ecog,2));
        audio_data((max_size+1):end) = [];
        ecog(:,(max_size+1):end) = [];
    else
        'ERRORR!!!'
    end
end
    

%% Clip ends of audio and ecog
%% ***Shift audio relative to ecog - clip start and end ecog [=
%% Assume data vectors are the same length (double check this assumption);
Symm_Shift_Range = 1;   % Symmetric Time(s) of shift before & after
clip_size = fs*Symm_Shift_Range;
ecog = ecog(:,(1+clip_size):(end-clip_size));
for j = 1:size(ecog,1) % Channel
    ecog_ch = ecog(j,:);
    for i = 1:(2*clip_size) % Time
        audio_range = audio_data(i:(i+length(audio_data)-2*clip_size-1));
        r(j,i) = corr(audio_range, ecog_ch'); %,'type', 'Spearman');
    end
end

lag = (1:(2*clip_size))/(fs) - Symm_Shift_Range;

figure; plot(repmat(lag',1,length(channels)),r')
hold
legend('Ch 41-49','Ch 33-41','Ch 25-33','Ch 17-25')
xlabel('Audio Offset from ECoG (s)')
ylabel('Pearson Correlation Coefficient')
axis('tight')
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0], get(gca,'ylim'),'k','LineWidth',1.5)
if use_HG
    title('Cross Correlation Between ECoG High Gamma and Audio Envelope Over Time')
else
    title('Cross Correlation Between Broadband ECoG and Audio Envelope Over Time')
end
a =1;
end
