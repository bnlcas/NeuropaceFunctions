function [] = RNS_ERP_analysis(root_dir)
%% This script will load and process ECoG and VNS data from EC71, RNS1 and RNS2
% and generate a series of plots.
% This will form ERPs for both high gamma and broadband signal for sentence
% and phoneme erps, and plot spectral differences between the various channels
% for a given electrode over time Spectral differences
% and finally, plot a heat map of intensity of trial vs time time for each electrode.
%
% This will require the following Steps:
%
%% Load ECoG TIMIT Data (EC 71)
% Create BiPolar paired ERPs for Sentence onsets in Broadband
% Create BiPolar paired ERPs for Sentence onsets in HighGamma
% Create BiPolar paired ERPs for Phonemen onsets in Broadband
% Create BiPolar paired ERPs for Phonemen onsets in HighGamma
% Plot BiPolar paired Spectra of Relevant Channesl in ECoG
% Plot Heatmap of relevant channels
%
%% Load RNS data from Fall (RNS_1)
% Time lock data for relevant bands.
% Create ERPs for Sentence onsets in Broadband
% Create ERPs for Sentence onsets in HighGamma
% Create ERPs for Phonemen onsets in Broadband
% Create ERPs for Phonemen onsets in HighGamma
% Plot spectra
% Plot Heatmap of erps
%
%% Load RNS data from Spring (RNS_2)
% Time lock data for relevant bands.
% Create ERPs for Sentence onsets in Broadband
% Create ERPs for Sentence onsets in HighGamma
% Create ERPs for Phonemen onsets in Broadband
% Create ERPs for Phonemen onsets in HighGamma
% Plot spectra
% Plot Heatmap of erps
%
%% OUTPUTs
% Plot of Sentence Onset ERPs for 3 blocks in Broadband
% Plot of Sentence onset ERPs for 3 block in High Gamma
% Plot of Phoneme Onset ERPs for 3 blocks in Broadband
% Plot of Phoneme onset ERPs for 3 block in High Gamma
% Plot of Full Spectra for channels in 3 blocks
% Plot of Heatmap ERPs for 3 Blocks
%
%% Input:
% root_dir - the directory containing the data and functions related to the
% neuropace recordings.
%
%% Ex: 
% RNS_ERP_analysis('/Users/changlab/Documents/NeuroPace_Analysis')

addpath(genpath(root_dir))

flatten = @(x) x(:);
%% Start:
plot_sentence_erps = false;
plot_sent_HG = true;
plot_sent_BB = true;

plot_phn_erps = false;
plot_phn_HG = true;
plot_phn_BB = true;

plot_sentence_heat_maps = false;
plot_hm_hg = true;
plot_hm_bb = true;

plot_spectra = true;
plot_spectra_seperate = false;
plot_spectra_combined = true;
plot_speech_spectral_compare = true;

plot_psi = false;
plot_psi_bb = false;
plot_psi_hg = false;
plot_psi_sequence = false;

plot_audio_env_corr_hg = true;
plot_audio_env_corr_bb = true;
plot_audio_env_corr_single_set = true; % plots all the audio envelope correlation on 1 graph

plot_spectral_erps = true;

auto_save_figs = true;
if auto_save_figs
    save_dir = [root_dir '/Graphics/RNS_Summary_Figures_1000'];
    mkdir(save_dir)
end



%% load evnts matricies (contain auto TIMIT time synch data)
evnt_dir = [root_dir '/data/RNS_evnt_mats/'];

load([evnt_dir 'EC71_rnsform.mat']); % Grid recording event
evnt_ec71 = evnt;

load([evnt_dir 'RNS_fall_evnt.mat']); % First RNS event (fall)
evnt_rns_1 = evnt;

load([evnt_dir 'RNS_spring_evnt.mat']); % Second RNS reconding event (spring)
evnt_rns_2 = evnt;

%% Get offsets_blocks
%load([evnt_dir 'offsets_1.mat']);
%load([evnt_dir 'offsets_2.mat']);

%% Clip non-data and extra blocks
relevant_evnts_ec71 = [1:248];
evnt_ec71 = evnt_ec71(relevant_evnts_ec71);

%relevant_evnts_rns1 = [1:100, 151:250, 301:499];
relevant_evnts_rns1 = [1:250, 301:499]; % Block 6 was clipped...
%relevant_evnts_rns1 = [1:499];
evnt_rns_1 = evnt_rns_1(relevant_evnts_rns1);

%relevant_evnts_rns2 = [1:699];
relevant_evnts_rns2 = [51:549];
evnt_rns_2 = evnt_rns_2(relevant_evnts_rns2);

%% Transform directories to local disc
[evnt_ec71] = transform_evnt_dir(evnt_ec71, root_dir); % Necessary to re route directories:
[evnt_rns_1] = transform_evnt_dir(evnt_rns_1, root_dir);
[evnt_rns_2] = transform_evnt_dir(evnt_rns_2, root_dir);

%% Make Sentence ERPs:
if plot_sentence_erps
    %% EC 71
%    [sent_bb_ec71, time_axis_ec71bb] = make_sentence_erps_grid(evnt_ec71, false);
%    [sent_hg_ec71, time_axis_ec71hg] = make_sentence_erps_grid(evnt_ec71, true);
    [sent_bb_ec71, time_axis_ec71bb] = make_sentence_erps_vns(evnt_ec71, false,400);
    [sent_hg_ec71, time_axis_ec71hg] = make_sentence_erps_vns(evnt_ec71, true,400);
  
    %% RNS 1
    [sent_hg_rns1, time_axis_rns] = make_sentence_erps_vns(evnt_rns_1, true,250);
    [sent_bb_rns1, time_axis_rns] = make_sentence_erps_vns(evnt_rns_1, false,250);
    
    %% RNS 2:  
    [sent_hg_rns2, time_axis_rns] = make_sentence_erps_vns(evnt_rns_2, true, 250);
    [sent_bb_rns2, time_axis_rns] = make_sentence_erps_vns(evnt_rns_2, false, 250);

end

%% Make Long timescale ERPS (-1, 2 seconds)
if plot_spectral_erps
    [sent_bb_ec71_long, time_axis_ec71bb_long] = make_sentence_erps_vns(evnt_ec71, false,400,[-1 2]);

    [sent_bb_rns1_long, time_axis_rns_long] = make_sentence_erps_vns(evnt_rns_1, false,250,[-1 2]);

    [sent_bb_rns2_long, time_axis_rns_long] = make_sentence_erps_vns(evnt_rns_2, false,250,[-1 2]);
end



%% Make Phoneme ERPs
if plot_phn_erps
    %% EC 71
%    phn_bb_ec71 = make_phoneme_erps_grid(evnt_ec71, false);
%    phn_hg_ec71 = make_phoneme_erps_grid(evnt_ec71, true);
% 
%   for x = 1:10
%       disp(x)
%   end
% 
    phn_bb_ec71 = make_phoneme_erps_rns(evnt_ec71, false ,400);
    phn_hg_ec71 = make_phoneme_erps_rns(evnt_ec71, true, 400);


    %% RNS 1
    phn_bb_rns1 = make_phoneme_erps_rns(evnt_rns_1, false, 250);
    phn_hg_rns1 = make_phoneme_erps_rns(evnt_rns_1, true, 250);
    
    %% RNS 2
    phn_bb_rns2 = make_phoneme_erps_rns(evnt_rns_2, false, 250);
    phn_hg_rns2 = make_phoneme_erps_rns(evnt_rns_2, true, 250);
end

%% if auto_save_figs
if auto_save_figs
    if plot_sent_HG
        sent_HG_dir = [save_dir '/Sentence ERP HighGamma/'];
        mkdir(sent_HG_dir)
    end
    if plot_sent_BB
        sent_BB_dir = [save_dir '/Sentence ERP Broadband/'];
        mkdir(sent_BB_dir)
    end
    if plot_phn_HG
        phn_HG_dir = [save_dir '/Phoneme ERPs HighGamma/'];
        mkdir(phn_HG_dir)
    end
    if plot_phn_BB
        phn_BB_dir = [save_dir '/Phoneme ERPs Broadband/'];
        mkdir(phn_BB_dir)
    end
    if plot_spectra
        spectra_dir = [save_dir '/Spectrograms/'];
        mkdir(spectra_dir)
    end
    if plot_audio_env_corr_hg | plot_audio_env_corr_bb
        audio_env_xcorr_dir = [save_dir '/Audio Envelope XCorr/'];
        mkdir(audio_env_xcorr_dir)
        if plot_audio_env_corr_bb
            audio_env_xcorr_bb_dir = [audio_env_xcorr_dir '/Broadband/'];
            mkdir(audio_env_xcorr_bb_dir)
        end
        if plot_audio_env_corr_hg
            audio_env_xcorr_hg_dir = [audio_env_xcorr_dir '/High Gamma/'];
            mkdir(audio_env_xcorr_hg_dir)
        end
    end
    if plot_spectral_erps
        spectral_erp_dir = [save_dir '/Spectral_ERPs/'];
        mkdir(spectral_erp_dir);
    end
end




%% Plots:
if plot_sentence_erps
    if plot_sent_HG
        %plot_rns_chans(sent_hg_ec71 , time_axis_ec71hg, true, 'EC 71 Sentence Onset', [-0.35 1.15]);
        multi_plot_rns_chans('EC 71 Sentence Onset', [-0.35 2], sent_hg_ec71 , time_axis_ec71hg);
        z1 = gcf;
        
        %plot_rns_chans(sent_hg_rns1, time_axis_rns, true,'RNS 1 Sentence', [-0.35 1.15]);
        multi_plot_rns_chans('RNS 1 Sentence Onset', [-0.35 1.0], sent_hg_rns1, time_axis_rns);
        z2 = gcf;

        %plot_rns_chans(sent_hg_rns2, time_axis_rns, true,'RNS 2 Sentence', [-0.35 1.15]);
        multi_plot_rns_chans('RNS 2 Sentence Onset', [-0.35 1.4], sent_hg_rns2, time_axis_rns);
        z3 = gcf;

        multi_plot_rns_chans('High Gamma Sentence ERPs', [-0.35 2], sent_hg_ec71 , time_axis_ec71hg, sent_hg_rns1, time_axis_rns, sent_hg_rns2, time_axis_rns)
        z4 = gcf;
        
        if auto_save_figs
            savefig(z1, [sent_HG_dir 'EC71_Sentence_Onset_HG.fig'])
            savefig(z2, [sent_HG_dir 'RNS_1_Sentence_Onset_HG.fig'])
            savefig(z3, [sent_HG_dir 'RNS_2_Sentence_Onset_HG.fig'])
            savefig(z4, [sent_HG_dir 'Combined_Sentence_Onset_HG.fig'])
            close(z1, z2, z3, z4)
        end
    end
    if plot_sent_BB
        %plot_rns_chans(sent_bb_ec71 , time_axis_ec71bb, false, 'EC 71 Sentence Onset', [-0.8 0.6]);
        multi_plot_rns_chans('EC71 Sentence Onset ERPs', [-0.5 0.45], sent_bb_ec71 , time_axis_ec71bb);
        z1 = gcf;
        
        %plot_rns_chans(sent_bb_rns1, time_axis_rns, false,'RNS 1 Sentence', [-0.4 0.3]);
        multi_plot_rns_chans('RNS 1 Sentence Onset ERPs', [-0.4 0.34], sent_bb_rns1, time_axis_rns);
        z2 = gcf;
        
        %plot_rns_chans(sent_bb_rns2, time_axis_rns, false,'RNS 2 Sentence', [-0.5 0.4]);
        multi_plot_rns_chans('RNS 2 Sentence Onset ERPs', [-0.5 0.4], sent_bb_rns2, time_axis_rns);
        z3 = gcf;
        
        multi_plot_rns_chans('Broadband Sentence ERPs', [-0.5 0.45], sent_bb_ec71 , time_axis_ec71bb, sent_bb_rns1, time_axis_rns, sent_bb_rns2, time_axis_rns)
        z4 = gcf;
        
        if auto_save_figs
            savefig(z1, [sent_BB_dir 'EC71_Sentence_Onset_BB.fig'])
            savefig(z2, [sent_BB_dir 'RNS_1_Sentence_Onset_BB.fig'])
            savefig(z3, [sent_BB_dir 'RNS_2_Sentence_Onset_BB.fig'])
            savefig(z4, [sent_BB_dir 'Combined_Sentence_Onset_BB.fig'])
            close(z1, z2, z3, z4)
        end
    
    end
end

if plot_phn_erps
    fricatives = {'f', 'v', 'dh', 'th', 'g', 's', 'z', 'sh', 'ch', 'jh'};    
    plosives = {'b', 'g', 'p', 't', 'k', 'd'};
    nasals = {'n', 'm', 'ng'};
    vowels = {'aa', 'ao', 'ow', 'ah', 'ax', 'uh', 'ux', 'uw', 'iy', 'ih', ...
        'ix', 'ey', 'eh', 'ae', 'aw', 'ay'};

    if plot_phn_HG
        Plot_Phoneme_ERPs(phn_hg_ec71, true, 'EC 71', [-0.05 1.85], is_consonant_type(phn_hg_ec71, fricatives) , is_consonant_type(phn_hg_ec71, plosives), is_consonant_type(phn_hg_ec71, nasals), is_consonant_type(phn_hg_ec71, vowels));
        z1 = gcf;
        
        Plot_Phoneme_ERPs(phn_hg_rns1, true, 'RNS 1', [-0.5 0.55], is_consonant_type(phn_hg_rns1, fricatives) , is_consonant_type(phn_hg_rns1, plosives), is_consonant_type(phn_hg_rns1, nasals), is_consonant_type(phn_hg_rns1, vowels));
        z2 = gcf;
        
        Plot_Phoneme_ERPs(phn_hg_rns2, true, 'RNS 2', [-0.15 1.1], is_consonant_type(phn_hg_rns2, fricatives) , is_consonant_type(phn_hg_rns2, plosives), is_consonant_type(phn_hg_rns2, nasals), is_consonant_type(phn_hg_rns2, vowels));        
        z3 = gcf;
        
        if auto_save_figs
            savefig(z1, [phn_HG_dir 'EC71_Phoneme_ERPs_HG.fig'])
            savefig(z2, [phn_HG_dir 'RNS_1_Phoneme_ERPs_HG.fig'])
            savefig(z3, [phn_HG_dir 'RNS_2_Phoneme_ERPs_HG.fig'])
            close(z1, z2, z3)
        end
    
    end
    if plot_phn_BB
        Plot_Phoneme_ERPs(phn_bb_ec71, false, 'EC 71', [-.2 0.15], is_consonant_type(phn_bb_ec71, fricatives) , is_consonant_type(phn_bb_ec71, plosives), is_consonant_type(phn_bb_ec71, nasals), is_consonant_type(phn_bb_ec71, vowels));
        z1 = gcf;
        
        Plot_Phoneme_ERPs(phn_bb_rns1, false, 'RNS 1', [-0.2 0.16], is_consonant_type(phn_bb_rns1, fricatives) , is_consonant_type(phn_bb_rns1, plosives), is_consonant_type(phn_bb_rns1, nasals), is_consonant_type(phn_bb_rns1, vowels));
        z2 = gcf;
        
        Plot_Phoneme_ERPs(phn_bb_rns2, false, 'RNS 2', [-0.2 0.16], is_consonant_type(phn_bb_rns2, fricatives) , is_consonant_type(phn_bb_rns2, plosives), is_consonant_type(phn_bb_rns2, nasals), is_consonant_type(phn_bb_rns2, vowels));
        z3 = gcf;    
        
        if auto_save_figs
            savefig(z1, [phn_BB_dir 'EC71_Phoneme_ERPs_BB.fig'])
            savefig(z2, [phn_BB_dir 'RNS_1_Phoneme_ERPs_BB.fig'])
            savefig(z3, [phn_BB_dir 'RNS_2_Phoneme_ERPs_BB.fig'])
            close(z1, z2, z3)
        end
    end
end

if plot_sentence_heat_maps
    if plot_hm_hg
        %% High gamma Heat Maps
            %% EC 71
            figure; hold on;
            for i = 1:4
                subplot(1,4,i)
                plot_erp_heatmap(squeeze(sent_hg_ec71(i,:,:)), time_axis_ec71hg, [-3, 8], evnt_ec71);
                title(['EC 71 High Gamma ERP Heatmap in Ch ' num2str(i)])           
            end

            %% RNS 1
            figure; hold on;
            for i = 1:4
                subplot(1,4,i);
                plot_erp_heatmap(squeeze(sent_hg_rns1(i,:,:)), time_axis_rns, [-2 5]); %, evnt_rns_1);
                title(['RNS 1 High Gamma ERP Heatmap in Ch ' num2str(i)])
            end
            
            %% RNS 2
            figure; hold on;
            for i = 1:4
                subplot(1,4,i);
                plot_erp_heatmap(squeeze(sent_hg_rns2(i,:,:)), time_axis_rns, [-2 5]); %,evnt_rns_2);
                title(['RNS 2 High Gamma ERP Heatmap in Ch ' num2str(i)])
            end            
            
    end
    if plot_hm_bb
        %% Plot BROADBAND Heatmap ERPs
            %% EC 71
            figure; hold on;
            for i = 1:4
                subplot(1,4,i)
                plot_erp_heatmap(squeeze(sent_bb_ec71(i,:,:)), time_axis_ec71bb, [-4 4]); %, evnt_ec71);
                title(['EC 71 Broadband ERP Heatmap in Ch ' num2str(i)])           
            end

            %% RNS 1
            figure; hold on;
            for i = 1:4
                subplot(1,4,i);
                plot_erp_heatmap(squeeze(sent_bb_rns1(i,:,:)), time_axis_rns,[-0.8 0.8]); %, evnt_rns_1);
                title(['RNS 1 Broadband ERP Heatmap in Ch ' num2str(i)])
            end
            
            %% RNS 2
            figure; hold on;
            for i = 1:4
                subplot(1,4,i);
                plot_erp_heatmap(squeeze(sent_bb_rns2(i,:,:)), time_axis_rns,[-0.8 0.8]); %, evnt_rns_2);
                title(['RNS 2 Broadband ERP Heatmap in Ch ' num2str(i)])
            end            
            
    end
end

%% PLOT Spectrograms
if plot_spectra
    %% Plot Separate Spectra for each Event
    if plot_spectra_seperate
        plot_spectrogram_rns_evnt(evnt_ec71, 'EC 71', 300, [20 3000])

        plot_spectrogram_rns_evnt(evnt_rns_1, 'RNS 1', 500, [80 2000])

        plot_spectrogram_rns_evnt(evnt_rns_2, 'RNS 2', 500, [90 2000])
    end
    %% Plot every Event Spectrum on a single set of plots:
    if plot_spectra_combined
        plot_all = true;
        plot_rns = true;
<<<<<<< HEAD
        normalize_spectra = false; % Flag to Normalize the spectrograms
        if plot_all
            plot_spectrogram_rns_evnt(evnt_ec71, 'EC 71', 800, [], normalize_spectra, false)

            plot_spectrogram_rns_evnt(evnt_rns_1, 'RNS 1', 800, [], normalize_spectra, true)

            plot_spectrogram_rns_evnt(evnt_rns_2, 'RNS 2', 800, [], normalize_spectra, true)
=======
        if plot_all
            plot_spectrogram_rns_evnt(evnt_ec71, 'EC 71', 800, [], true, false)

            plot_spectrogram_rns_evnt(evnt_rns_1, 'RNS 1', 800, [], true, true)

            plot_spectrogram_rns_evnt(evnt_rns_2, 'RNS 2', 800, [], true, true)
>>>>>>> e96c2f3a7b4a9c501e2ae4c343d4b60513720af8
            for i = 1:4
                subplot(1,4,i);
                ylim([10^-3 3*10^1])
                legend('EC 71', 'RNS 1', 'RNS 2')
                title('Power Normalized Spectral Density')
            end
            z1 = gcf;
        end
        if plot_rns
<<<<<<< HEAD
            plot_spectrogram_rns_evnt(evnt_rns_1, 'RNS 1', 800, [], normalize_spectra, false)

            plot_spectrogram_rns_evnt(evnt_rns_2, 'RNS 2', 800, [], normalize_spectra, true)
            for i = 1:4
                subplot(1,4,i);
                if normalize_spectra
                    ylim([10^-2 1.2*10^1])
                    ylabel('Normalized Power Density (P^2/Hz)');
                    title('Normalized Spectral Power Density')
                else
                    ylim([10^-1 8*10^1])
                    ylabel('Power Density (P^2/Hz)');
                    title('Spectral Power Density')
                end
                legend('RNS 1', 'RNS 2')               
=======
            plot_spectrogram_rns_evnt(evnt_rns_1, 'RNS 1', 800, [], true, false)

            plot_spectrogram_rns_evnt(evnt_rns_2, 'RNS 2', 800, [], true, true)
            for i = 1:4
                subplot(1,4,i);
                ylim([10^-2 1.2*10^1])
                legend('RNS 1', 'RNS 2')
                title('Power Normalized Spectral Density')
>>>>>>> e96c2f3a7b4a9c501e2ae4c343d4b60513720af8
            end
            z2 = gcf;
        end
        if plot_speech_spectral_compare
            [ecog_speech, ecog_nonspeech] = gather_raw_speech_nonspeech_ecog(evnt_rns_1);
            plot_spectrum_chans(ecog_speech, 250, 'Power Normalized Spectrum of Speech/Non-Speech in RNS 1', 800, [], true, false);
            plot_spectrum_chans(ecog_nonspeech, 250, 'Power Normalized Spectrum of Speech/Non-Speech in RNS 1', 800, [], true, true);
            for i = 1:4
                subplot(1,4,i)
                ylim([0.01 15])
                legend('Speech', 'Non-Speech')
                title({'Power Normalized Spectral Density of'; 'Speech/Non-Speech in RNS 1'})
            end
            z3 = gcf;
            
            [ecog_speech, ecog_nonspeech] = gather_raw_speech_nonspeech_ecog(evnt_rns_2);
            plot_spectrum_chans(ecog_speech, 250, 'Power Normalized Spectrum of Speech/Non-Speech in RNS 2', 800, [], true, false);
            plot_spectrum_chans(ecog_nonspeech, 250, 'Power Normalized Spectrum of Speech/Non-Speech in RNS 2', 800, [], true, true);
            for i = 1:4
                subplot(1,4,i)
                legend('Speech', 'Non-Speech')
                ylim([0.01 11])
                title({'Power Normalized Spectral Density of'; 'Speech/Non-Speech in RNS 2'})
            end
            z4 = gcf;
            
            [ecog_speech, ecog_nonspeech] = gather_raw_speech_nonspeech_ecog(evnt_ec71);
            plot_spectrum_chans(ecog_speech, 400, 'Power Normalized Spectrum of Speech/Non-Speech in EC71', 700, [], true, false);
            plot_spectrum_chans(ecog_nonspeech, 400, 'Power Normalized Spectrum of Speech/Non-Speech in EC71', 700, [], true, true);
            for i = 1:4
                subplot(1,4,i)
                legend('Speech', 'Non-Speech')
                ylim([10^-3 3*10^1])
                title({'Power Normalized Spectral Density of'; 'Speech/Non-Speech in EC 71'})
            end
            z5 = gcf;
        end
        
        if auto_save_figs
            savefig(z1, [spectra_dir 'EC 71 & RNS 1 & RNS 2 spectra.fig'])
            savefig(z2, [spectra_dir 'RNS 1 & RNS 2 spectra.fig'])
            close(z1, z2)
            if plot_speech_spectral_compare
                savefig(z3, [spectra_dir, 'Comparison of RNS 1 spectra with without Speech.fig'])
                savefig(z4, [spectra_dir, 'Comparison of RNS 2 spectra with without Speech.fig'])
                savefig(z5, [spectra_dir, 'Comparison of EC 71 spectra with without Speech.fig'])

                close(z3, z4, z5)
            end
        end

    end
end

%% Plot PSI's
if plot_psi & plot_phn_erps
    %% Plot High Gamma PSI's
    if plot_psi_hg
        tmp = PSI_RNS(phn_hg_ec71 ,400, false, 0.15, 0.1);
        title('PSI for EC 71 EMU Recording (High Gamma)');
        caxis([0 10])

        tmp = PSI_RNS(phn_hg_rns1,250,true, 0.15, 0.1);
        title('PSI for RNS Recording 1 (High Gamma)');
        caxis([0 10])

        tmp = PSI_RNS(phn_hg_rns2,250,true, 0.15, 0.1);
        title('PSI for RNS Recording 2 (High Gamma)')
        caxis([0 10])
    end
    %% Plot BB PSI's
    if plot_psi_bb
        tmp = PSI_RNS(phn_bb_ec71 ,400, true, 0.15, 0.1);
        title('PSI for EC 71 EMU Recording (Broadband)');
        caxis([0 20])

        tmp = PSI_RNS(phn_bb_rns1,250, true, 0.15, 0.1);
        title('PSI for RNS Recording 1 (Broadband)');
        caxis([0 20])

        tmp = PSI_RNS(phn_bb_rns2,250,true, 0.15, 0.1);
        title('PSI for RNS Recording 2 (Broadband)')
        caxis([0 20])
    end
    if plot_psi_sequence
        psi_tc_ec71_bb = psi_timecourse(phn_bb_ec71, 'BB EC 71');
      
        psi_tc_rns1_bb = psi_timecourse(phn_bb_rns1, 'BB RNS 1');
        
        psi_tc_rns2_bb = psi_timecourse(phn_bb_rns2, 'BB RNS 2');
        
        
        psi_tc_ec71_hg = psi_timecourse(phn_hg_ec71, 'HG EC 71');
      
        psi_tc_rns1_hg = psi_timecourse(phn_hg_rns1, 'HG RNS 1');
        
        psi_tc_rns2_hg = psi_timecourse(phn_hg_rns2, 'HG RNS 2');
        
    end

end

%% Plot the Correlation between the amplitude of the audio and ECoG
if plot_audio_env_corr_hg
    %x_corr_ec71 = audio_env_ecog_corr(evnt_ec71, true, false);
    x_corr_waves_ec71 = audio_env_ecog_corr_grid(evnt_ec71, true,false,[-0.05 0.08],'EC 71');
    z1 = gcf;
    
    x_corr_rns1 = audio_env_ecog_corr(evnt_rns_1, true, true, [-0.02, 0.15], 'RNS 1');
    z2 = gcf;
    
    x_corr_rns2 = audio_env_ecog_corr(evnt_rns_2, true, true, [-0.05, 0.23], 'RNS 2');
    z3 = gcf;
    
    if plot_audio_env_corr_single_set
        Plot_audio_env_corr_unified(x_corr_waves_ec71, x_corr_rns1, x_corr_rns2, 'TimeSpan',1, 'YScale', [-0.05 0.23]);
        z4 = gcf;
    end
    
    
    if auto_save_figs
        savefig(z1, [audio_env_xcorr_hg_dir 'EC 71 Audio-ECoG xcorr HG.fig'])
        savefig(z2, [audio_env_xcorr_hg_dir 'RNS 1 Audio-ECoG xcorr HG.fig'])
        savefig(z3, [audio_env_xcorr_hg_dir 'RNS 2 Audio-ECoG xcorr HG.fig'])
        close(z1, z2, z3)
        if plot_audio_env_corr_single_set
            savefig(z4, [audio_env_xcorr_hg_dir 'Combined Audio-ECoG xcorr HG.fig'])
            close(z4)
        end
    end

end
if plot_audio_env_corr_bb
     x_corr_waves_ec71_bb = audio_env_ecog_corr_grid(evnt_ec71, false ,false,[-0.05 0.05],'EC 71 Broadband'); 
     z1 = gcf;
     
     x_corr_rns1_bb = audio_env_ecog_corr(evnt_rns_1, false, true, [-0.1, 0.05], 'RNS 1 Broadband');
     z2 = gcf;
     
     x_corr_rns2_bb = audio_env_ecog_corr(evnt_rns_2, false, true, [-0.06, 0.03], 'RNS 2 Broadband');
     z3 = gcf;
     if plot_audio_env_corr_single_set
         Plot_audio_env_corr_unified(x_corr_waves_ec71_bb, x_corr_rns1_bb, x_corr_rns2_bb, 'TimeSpan',1, 'YScale', [-0.1 0.05]);
        z4 = gcf;
     end
     
     
     if auto_save_figs
        savefig(z1, [audio_env_xcorr_bb_dir 'EC 71 Audio-ECoG xcorr BB.fig'])
        savefig(z2, [audio_env_xcorr_bb_dir 'RNS 1 Audio-ECoG xcorr BB.fig'])
        savefig(z3, [audio_env_xcorr_bb_dir 'RNS 2 Audio-ECoG xcorr BB.fig'])
        close(z1, z2, z3)
        if plot_audio_env_corr_single_set
            savefig(z4, [audio_env_xcorr_bb_dir 'Combined Audio-ECoG xcorr BB.fig'])
            close(z4)
        end
    end


% old xlim for non-bb envelope
% % % % %     x_corr_waves_ec71_bb = audio_env_ecog_corr_grid(evnt_ec71, false ,false,[-0.05 0.04],'EC 71 Broadband'); 
% % % % %     x_corr_rns1_bb = audio_env_ecog_corr(evnt_rns_1, false, true, [-0.05, 0.04], 'RNS 1 Broadband');
% % % % %     x_corr_rns2_bb = audio_env_ecog_corr(evnt_rns_2, false, true, [-0.05, 0.04], 'RNS 2 Broadband');
end



if plot_spectral_erps
    plot_event_spectrograms_iii(sent_bb_ec71_long, time_axis_ec71bb_long);
    z1 = gcf;
    
    plot_event_spectrograms_iii(sent_bb_rns1_long, time_axis_rns_long);
    z2 = gcf;
    
    plot_event_spectrograms_iii(sent_bb_rns2_long, time_axis_rns_long);    z3 = gcf;
    z3 = gcf;
    
    if auto_save_figs
        savefig(z1, [spectral_erp_dir, 'EC 71 Spectral ERPs.fig'])
        savefig(z2, [spectral_erp_dir, 'RNS 1 Spectral ERPs.fig'])
        savefig(z3, [spectral_erp_dir, 'RNS 2 Spectral ERPs.fig'])
        close(z1, z2, z3)
    end
end




end

