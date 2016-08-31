function [] = Preprocess_EC71_grid_data(root_dir)
%This function creates loads the raw data for EC71 and converts it into 4
%channels of bipolar montaged RAW data and high gamma.
% THis is then saved in a specificed directory

addpath(genpath(root_dir));
%% INPUT:
% root_dir - is given as the directory containing NeuroPace_Analysis:
% EX: root_dir = '/Users/changlab/Documents/NeuroPace_Analysis'
load([root_dir '/data/RNS_evnt_mats/EC71.mat']);
root_dirs = {[root_dir '/data/EC71/EC71_B1'],...
    [root_dir '/data/EC71/EC71_B5']};

%% Define BiPolar Pairings
%channels = [17 25 33 41]; % Channels of interest
%ref_chans = [18 26 34 42]; % ORIGINAL SET


%channels = [42 18 41 17];
%ref_chans = [34 26 33 25]; %CORRECT BP polar?
channels = [42 26 41 25];


bipolar_ref = false;

%fs_ecog = 250; % Resample Everything to 250 Hz?
fs_ecog = 400;
apply_4hz_hp = true; % apply mimic high pass filter:
if apply_4hz_hp
          hpFilt = designfilt('highpassiir','FilterOrder',2, ...
         'PassbandFrequency',4,'PassbandRipple',0.1, ...
         'SampleRate',fs_ecog);
end


%% Loop through the relevant EC71 directories:
%erps_total = []
for dir = root_dirs
    num_chans = 64;
    rns_style_dir = [dir{1} '/RNS_form'];
    mkdir(rns_style_dir)
    
    tmp = strsplit(rns_style_dir,'/');
    tmp{end-1} = [tmp{end-1} '_RAW'];
    rns_style_dir_raw = strjoin(tmp,'/');
    mkdir(rns_style_dir_raw);
    
    %% load RAW Data:
    raw_path = [dir{1} '/RawHTK'];

    ecog_grid = [];
    for j = 1:num_chans
        ch_path = [raw_path, '/Wav2', num2str(j), '.htk'];

        [ch_data, fs_data] = readhtk(ch_path);
        ch_data_ds = resample(ch_data, 2^11, 5^6);
        %ch_data_ds = resample(ch_data,(5*2^11), (8*5^6));
        
% % %         resample(ecog.data(i,:),m*2^11, n*5^6);
% % %         (5*2^11)/(8*5^6)
% % %         
        if apply_4hz_hp
            ch_data_filt = filtfilt(hpFilt,ch_data_ds);
            ch_data_ds = ch_data_filt;
        end
        ch_data_ds_z = (ch_data_ds - mean(ch_data_ds))/std(ch_data_ds); % Zscore
        ecog_grid = [ecog_grid; ch_data_ds_z];
    end
    if mod(size(ecog_grid,2),2) == 1
        ecog_grid(:,end) = [];
    end
%% Take Common median reference - Adds noise to high gamma
% % % %     for j = 1:size(ecog_grid,2)
% % % %         ecog_grid(:,j) = ecog_grid(:,j) - median(ecog_grid(:,j));
% % % %     end
    %% Notch Filter 60 & 120 Hz:
    temp_struct.data = ecog_grid;
    temp_struct_out= applyLineNoiseNotch_60HzHarmonics(temp_struct, fs_ecog);
    ecog_grid = temp_struct_out.data;
    %ecog_grid = gdivide(gsubtract(ecog_grid, mean(ecog_grid,2)),std(ecog_grid,[],2));
% % % % % % %     erps = [];
% % % % % % %     TIMIT_block_delay = 0.12;
% % % % % % %     for j = 1:length(evnt)
% % % % % % %         if strcmpi(evnt(j).dpath, dir)
% % % % % % %             start_pt = round((evnt(j).StartTime+TIMIT_block_delay)*fs_ecog);
% % % % % % %             data = ecog_grid(:,(start_pt-round(fs_ecog/2)):(start_pt+round(fs_ecog)));
% % % % % % %             erps = cat(3, erps, data);
% % % % % % %         end
% % % % % % %     end
% % % % % % %     erps_total = cat(3,erps_total, erps);
% % % % % % %     time_axis = linspace(-0.5, 1, size(erps,2));
% % % % % % %                 
    
    if bipolar_ref
        ecog_bp_raw = ecog_grid(channels,:) - ecog_grid(ref_chans,:);
    else
 %       tmp = ecog_grid(ref_chans,:);
 %       tmp2 = ecog_grid(channels,:);
 %       channels = [ref_chans, channels];
 %       ecog_bp_raw = [tmp;tmp2];
        ecog_bp_raw = ecog_grid(channels,:);
    end
    ecog_bp_raw = gdivide(gsubtract(ecog_bp_raw, mean(ecog_bp_raw,2)),std(ecog_bp_raw,[],2));
    ecog_bp_hg = zeros(size(ecog_bp_raw));
    for j = 1:length(channels)
        hg_tmp = Hilbert_HG_std_ii(ecog_bp_raw(j,:));
        ecog_bp_hg(j,:) = (hg_tmp - mean(hg_tmp))/std(hg_tmp);
    end



    %% SAVE RNS Style BiPolar Montaged Data
    rns_style_dir_raw = [rns_style_dir_raw, '/HilbAA_70to150_8band'];
    mkdir(rns_style_dir_raw)
    rns_style_dir_Hilb = [rns_style_dir, '/HilbAA_70to150_8band'];
    mkdir(rns_style_dir_Hilb)

    for j = 1:length(channels)
        writehtk([rns_style_dir_raw, '/Wav1' num2str(j),'.htk'], ecog_bp_raw(j,:), fs_ecog);

        writehtk([rns_style_dir_Hilb, '/Wav1' num2str(j),'.htk'], ecog_bp_hg(j,:),fs_ecog);
    end
        
end


% % erps = erps_total;
% %     %Plot_Sub_Grid(time_axis, false, erps)
% %     grd = flipud(reshape(1:64,8,8))';
% %     bounds(2) = max(max(mean(erps,3) + nansem(erps,3)));
% %     bounds(1) = min(min(mean(erps,3) - nansem(erps,3)));
% %     figure;
% %     for j = 1:size(erps,1)
% %         p = plotGridPosition(j, 64, 8); 
% %         subplot('position',p);
% %         shadedErrorBar(time_axis, mean(squeeze(erps(grd(j),:,:)),2), nansem(squeeze(erps(grd(j),:,:)),2));
% %         hold on
% %         ylim(bounds)
% %         xlim([-0.2 1])
% %         plot([0 0], get(gca, 'YLim'),'k')
% %         axis tight;
% %         plot(get(gca, 'XLim'), [0 0], 'k')
% %         text(0,(max(get(gca,'YLim')) - (max(get(gca,'YLim')) * 0.1)),num2str(grd(j)))
% %     end

        
end
