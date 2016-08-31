%% Preprocess RNS files (does not generate plots, etc)

function [timit_num, is_clipped] = Preprocess_RNS_gen_ii(appelation, dat_dir, TIMIT_Blocknum, badfile)
% Takes raw data from 'dat_dir' and saves the preprocess, and assigns output
% folders with the tag designated by the inpiut 'appelation'.
% data is saved offset adjusted
% is_dropped is a boolean that is true for blocks that experienced clipped

% appelation = 'spring_2016'; % or 'fall_2015', or '0012',,,
% dat_dir = '/Users/changlab/Documents/changrepo/matlab/analysis/RNS/RNS Follow Up Spring 2016/UCSF_SM_3411999_2016-03-30_2016-03-30 EXTERNAL #PHI-selected';

% Num Sure About how to add TIMIT Blocknum

%% SPRING
% % % % % % %[offsets, timit_num, is_clipped] = Preprocess_RNS_gen_ii('spring', '/Users/changlab/Documents/changrepo/matlab/analysis/RNS/RNS Follow Up Spring 2016/UCSF_SM_3411999_2016-03-30_2016-03-30 EXTERNAL #PHI-selected');
% % % % % % TIMIT_Blocknum = [1:10, 1, 1, 1, 1]; % Totally arbitrary - based on physical event.
% % % % % % badfile = '131038325241940000.dat';

%[timit_num, is_clipped] = Preprocess_RNS_gen_ii('spring', '/Users/changlab/Documents/changrepo/matlab/analysis/RNS/RNS Follow Up Spring 2016/UCSF_SM_3411999_2016-03-30_2016-03-30 EXTERNAL #PHI-selected', [1:10, 1, 1, 1, 1], '131038325241940000.dat');



% % % %% FALL:
% % % [offsets, timit_num, is_clipped] = Preprocess_RNS_gen_ii('fall', '/Users/changlab/Documents/changrepo/matlab/analysis/RNS/RNS TIMIT Exp Data Files 111214/14464527_2014-11-12_21.14.19');
% % % TIMIT_Blocknum = [1:10];
% % % badfile = '';

%[timit_num, is_clipped] = Preprocess_RNS_gen_ii('fall', '/Users/changlab/Documents/changrepo/matlab/analysis/RNS/RNS TIMIT Exp Data Files 111214/14464527_2014-11-12_21.14.19', [1:10], '');


%% Establish Directories:
tmp = strsplit(dat_dir, '/');
tmp{end} = 'ECoG_Raw';
raw_dir = strjoin(tmp, '/');
mkdir(raw_dir);

data_out_dir = ['/Users/changlab/Documents/data/RNS_' appelation];
data_out_dir_raw = ['/Users/changlab/Documents/data/RNS_' appelation '_RAW'];
cd(dat_dir)
mkdir(data_out_dir)
mkdir(data_out_dir_raw)


dat_files = dir;

for i = 1:length(dat_files)
    file_name = dat_files(i);
    if substrcmp(lower(file_name.name), '.lay')
        lay_filename = file_name.name; % This doesn't really matter, any .lay file will work to extract ECoG Data
    end
end

%% Select good dat files:
i = 1;

while i <= length(dat_files)
    file_name = dat_files(i);
    if ~substrcmp(lower(file_name.name), '.dat') | strcmp(file_name.name, badfile)
        dat_files(i) = [];
        i = i-1;
    end
    i = i + 1;
end

dat_filenames = {''};
for i = 1:length(dat_files)
    file_name = dat_files(i);
    dat_filenames{i} = file_name.name;
end


%% Generate ERPs:

TIMIT_audio_dir = '/Users/changlab/Documents/changrepo/matlab/analysis/RNS/TIMIT_RNS';
cd(TIMIT_audio_dir)

is_clipped = false(length(TIMIT_Blocknum),1);
zero_clipped_dat = true; % Parameter to zero tracts of data that were dropped
clipped_val = -512;
% % % % % % %ecog_total = [];
% % % % % % %ecog_total_raw = [];
% % % % % % %ecog_total_std = [];
for i = 1:length(TIMIT_Blocknum)
    [ECoG_hdr, ECoG_data] = RNS_ReadECoGData(dat_filenames{i}, lay_filename);
    % Clean dropped data:
    is_clipped_pt = (ECoG_data(1,:) == clipped_val);
    is_clipped(i) = (sum(is_clipped_pt) > 0);
    if is_clipped(i) & zero_clipped_dat
         ECoG_data(:,is_clipped_pt) = repmat(median(ECoG_data(:,~is_clipped_pt),2),1, sum(is_clipped_pt));
    end
    fs_ecog = 250;
    
    % Remove 64 sample dropping artifact
    ECoG_data(1,:) = remove_64_samp_artifact(ECoG_data(1,:));
    ECoG_data(2,:) = remove_64_samp_artifact(ECoG_data(2,:));
    
    % notch filter data:
    %temp_struct.data = ECoG_data;
    %temp_struct_out= applyLineNoiseNotch_60HzHarmonics(temp_struct, fs_ecog);
    %ECoG_data = temp_struct_out.data;

    % Z Normalize ECoG:
    ECoG_data = gdivide(gsubtract(ECoG_data, mean(ECoG_data,2)),std(ECoG_data,[],2));
    
    
    % Amplify onset
    ecog_HG_AA = Hilbert_HG_transmax(ECoG_data);
    
    % RAW Data
    ecog_raw = ECoG_data;
    
    % High Gamma
    ecog_std = Hilbert_HG_std_ii(ECoG_data);
    ecog_std = gdivide(gsubtract(ecog_std, mean(ecog_std,2)),std(ecog_std,[],2));
    
    
    ECoG_data = ecog_HG_AA;
    [audio,fs] = audioread(['timit_block' num2str(TIMIT_Blocknum(i)) '.wav']);
    
    
    %% Block Audio into sentence Onsets:
    %audio_ds = decimate(audio,round(fs/fs_ecog));
    audio_diff = abs(diff(audio));
    thresh = 10^-5;
    is_on = audio_diff > thresh;
    is_onset = false(size(audio_diff));
    off_time = 500; % Number of points preciding  offset
    for j = (1+off_time):(length(is_onset)-1);
        if (sum(is_on((j-off_time):j)) == 0) & is_on(j+1)
            is_onset(j+1) = true;
        end
    end
    onset_inds = round(find(is_onset)*(fs_ecog/fs));
    
    %% Section ECoG Data by onset
    erp_dur = 10; % seconds on length for an ERP
    center_erp = true; % determine wether erp starts before or after audio
    center_offset = 2; % Shift from center
    ecog_erps = zeros(size(ECoG_data,1),erp_dur*fs_ecog,length(onset_inds));
    ecog_raw_erps = zeros(size(ECoG_data,1),erp_dur*fs_ecog,length(onset_inds));
    ecog_std_erps = zeros(size(ECoG_data,1),erp_dur*fs_ecog,length(onset_inds));
    
    null_trial = false(length(onset_inds),1);
    max_ind = size(ECoG_data,2);
    for j = 1:length(onset_inds)
        if center_erp
            start_ind = onset_inds(j)-round(0.5*erp_dur*fs_ecog - center_offset*fs_ecog);
            end_ind = start_ind+round(erp_dur*fs_ecog)-1;
            time_axis = center_offset+linspace(-0.5*erp_dur,0.5*erp_dur, erp_dur*fs_ecog);
        else
            start_ind = onset_inds(j);
            end_ind = onset_inds(j)+erp_dur*fs_ecog-1;
            time_axis = linspace(0, erp_dur, erp_dur*fs_ecog);
        end
        if end_ind < max_ind & start_ind > 0
             ecog_erps(:,:,j) = ECoG_data(:,start_ind:end_ind);
             ecog_raw_erps(:,:,j) = ecog_raw(:,start_ind:end_ind);
             ecog_std_erps(:,:,j) = ecog_std(:,start_ind:end_ind);
        else
            null_trial(j) = true;
        end
    end
    ecog_erps(:,:,null_trial) = [];
    ecog_raw_erps(:,:,null_trial) = [];
    ecog_std_erps(:,:,null_trial) = [];
    
    time_axis_diff = time_axis(1:(end-1)) + time_axis(find(time_axis>0,1))/2;

    plot_erps = false;
    if plot_erps
         figure;
         for j = 1:size(ecog_erps,1)
             subplot(1, size(ecog_erps,1), j)
             ecog_ch = squeeze(ecog_erps(j,:,:));
             shadedErrorBar(time_axis , mean(ecog_ch,2), nansem(ecog_ch,2))
             %ecog_ch = diff(squeeze(ecog_erps(j,:,:)));
             %shadedErrorBar(time_axis_diff , mean(ecog_ch,2), nansem(ecog_ch,2))
             axis tight
         end
    end
    
    sig_ch = 2;
    sig_ch_diff = diff(squeeze(ecog_erps(sig_ch,:,:)));
    [~,max_ind] = max(mean(diff(squeeze(ecog_erps(sig_ch,:,:))),2));
    %% SHIFT DATA back
    onset_erp_lag = 0.2; % time of delay of max derivate from stim onset (taken from EC71)
    max_ind = max_ind - round(onset_erp_lag*fs_ecog);
    
    shift_time = time_axis_diff(max_ind); % necessary to adjust for window shifint
    
    %% Clear Data before max in worth of space
    ecog_raw(:,1:round(shift_time*fs_ecog)) = [];
    ecog_std(:,1:round(shift_time*fs_ecog)) = [];
    
    
    
    %% Save Data from block
    
    data_out_dir_block = [data_out_dir ,'/B' num2str(i) '/HilbAA_70to150_8band/'];
    data_out_dir_block_raw = [data_out_dir_raw ,'/B' num2str(i) '/HilbAA_70to150_8band/'];
    data_out_dir_artifacts = [data_out_dir_raw, '/B' num2str(i) '/Artifacts/'];

    cd(data_out_dir_block)
    for j = 1:4
        writehtk([data_out_dir_block 'Wav1' num2str(j) '.htk'], ecog_std(j,:),fs_ecog);
        
        writehtk([data_out_dir_block_raw, 'Wav1' num2str(j) '.htk'], ecog_raw(j,:),fs_ecog);
    end
    
    %% Adjust bad_time segments
%     load([data_out_dir_artifacts 'badTimeSegments.mat'])
%     badTimeSegments = badTimeSegments - shift_time;
%     save([data_out_dir_artifacts 'badTimeSegments_shift.mat'], 'badTimeSegments')
end
    

timit_num = TIMIT_Blocknum;
end
