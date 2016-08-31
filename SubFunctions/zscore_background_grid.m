function [zscored_data] = get_background_zscore(input_data, use_HG)
%% Takes a n chann x time_pts x trials object and zsocres

% channels = [41 33 25 17]; % Channels of interest
% ref_chans = [49 41 33 25];
% relevant_channels = [channels, ref_chans];

%% Get Directories:
dpaths = {''};
for i = 1:length(evnt)
    dpaths{i} = evnt(i).dpath;
end
dpaths = unique(dpaths)

for i = 1:length(dpaths)
if use_HG
        data_path = [data_path, '/HilbAA_70to150_8band'];
    else
        data_path_split = strsplit(data_path, '/'); 
        data_path_split{end-1} = [data_path_split{end-1} '_RAW'];

        data_path = strjoin(data_path_split,'/');
        data_path = [data_path '/HilbAA_70to150_8band'];
    end


for i = 1:size(input_data,1)
    if use_HG
        data_path_B1 = '/Users/changlab/Documents/data/EC71/EC71_B1/HilbAA_70to150_8band_BP';
        data_path_B5 = '/Users/changlab/Documents/data/EC71/EC71_B5/HilbAA_70to150_8band_BP';
    else
        data_path_B1 = '/Users/changlab/Documents/data/EC71/EC71_B1/RawHTK';
        data_path_B5 = '/Users/changlab/Documents/data/EC71/EC71_B5/RawHTK';
    end
    
    if use_HG
        data1 = mean(readhtk([data_path_B1, '/Wav2', num2str(i),'.htk']),1);
        data5 = mean(readhtk([data_path_B5, '/Wav2', num2str(i),'.htk']),1);
        data = [data1, data5];
    else
        data1 = readhtk([data_path_B1, '/Wav2', num2str(i),'.htk']);
        %data1 = decimate(data1,12);
        data5 = readhtk([data_path_B5, '/Wav2', num2str(i),'.htk']);
        %data5 = decimate(data5,12);
        data = [data1, data5];
    end
    mean_ch(i) = mean(data);
    std_ch(i) = std(data);
end

for i = 1:size(input_data,1)
    ch_data = squeeze(input_data(i,:,:));
    zscored_data(i,:,:) = (ch_data-mean_ch(i))/std_ch(i);
end

end
