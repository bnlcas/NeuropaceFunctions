function ecog_HG_AA = Hilbert_HG_transmax(ecog)
%% Creates an output measure of the intensity of the High Gamma
%% From the HilberTransform
%% Attempts to Maximize the transient repsonse electrodes:

%% Calculates 
fs_ecog = 250;
freq_axis = (fs_ecog/2)*(0:floor(size(ecog,2)/2-1))/floor(size(ecog,2)/2-1);
freq_axis = [freq_axis, sort(freq_axis,'descend')];

freq_samples = [65.19885, 71.985375, 79.478298, 87.75115, 96.885131, 106.969857, 118.1042979];
window_widths = [3.1490865, 3.308923, 3.476873, 3.653347, 3.838779, 4.0336231, 4.23835625];
%window_widths = window_widths;


freq_samples = 55:0.5:120;
window_widths = 1.5*ones(size(freq_samples));


ecog_hg = zeros(size(ecog));
ecog_HG_AA = zeros([size(ecog),length(freq_samples)]);


for k = 1:length(freq_samples)

    %% Filter ECoG Channels
    gaussian_window = @(f) exp(-((f-freq_samples(k)).^2)/(2*window_widths(k).^2));
    for i = 1:size(ecog,1)
        ecog_ch = ecog(i,:);
        ecog_fft = fft(ecog_ch);
        
        T=length(ecog_ch);
        h = zeros(1,T);
        if 2*fix(T/2)==T %if T is even
            h([1 T/2+1]) = 1;
            h(2:T/2) = 2;
        else
            h(1) = 1; h(2:(T+1)/2) = 2;
        end

        gauss_filter = gaussian_window(freq_axis);
        
        hilbdata=ifft(ecog_fft(end,:).*(gauss_filter.*h),T);
        envData=abs(hilbdata);
 
        ecog_HG_AA(i,:,k) =envData;
 %       ecog_fft = gauss_filter.*ecog_fft;
  %      ecog_hg(i,:) = real(ifft(ecog_fft));

    end


    %     fourier_power_plot(mean(ecog_hg,1), fs_ecog);
%     for i = 1:size(ecog_hg,1);
%         ecog_HG_AA(i,:,k) = abs(hilbert(ecog_hg(i,:)));
%     end
    

end

%% Z-score each band and Reduce Dimension
% for i = 1:size(ecog_HG_AA,1)
%     for j = 1:size(ecog_HG_AA,3)
%         ecog_HG_AA(i,:,j) = (squeeze(ecog_HG_AA(i,:,j))-mean(squeeze(ecog_HG_AA(i,:,j))))/std(squeeze(ecog_HG_AA(i,:,j)));
%     end
% end

ecog_HG_AA = mean(ecog_HG_AA,3);




% 
% 
% for i = 1:size(ecog_HG_AA,1)
%     test_mean =  squeeze(mean(squeeze(ecog_HG_AA(i,:,:)),2))';
%     for j = 1:size(ecog_HG_AA,3)
%         ecog_HG_AA(i,:,j) = squeeze(ecog_HG_AA(i,:,j)) - test_mean;
%     end
% end
% 
% 
% 
% tmp1 = zeros(size(ecog_HG_AA,1),size(ecog_HG_AA,2));
% for i = 1:size(ecog_HG_AA,1)
%     pca_mat = pca(squeeze(ecog_HG_AA(i,:,:)));
%     tmp1(i,:)= squeeze(ecog_HG_AA(i,:,:))*pca_mat(:,1);
% end
% ecog_HG_AA = tmp1;
% 
% tmp1 = zeros(size(ecog_HG_AA,1),size(ecog_HG_AA,2));
% tmp2= zeros(size(ecog_HG_AA,2),size(ecog_HG_AA,3));
% % Take PCA of resulting AA for each channel
% for i = 1:size(ecog_HG_AA,1)
%     tmp2 = squeeze(ecog_HG_AA(i,:,:));
%     pca_mat = pca(tmp2);
%     tmp1(i,:) = tmp2*pca_mat(:,1);
% end
% ecog_HG_AA = tmp1;
% 
% 
% 
% 
