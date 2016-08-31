function smap = PSI_GRID(Phn, fs)
% function smap = PSI(dpath, subject)
%
% Get the PSI (phoneme selectivity index) matrix (smap) for a given subject
% 
% Inputs:
%       dpath [string]: path to your data file
%       subject [str]: your subject (e.g. 'EC56')
%       fs [int]: sampling rate of the data
%
% Need a file e.g. D_20bef50after_start.mat
% with variables
%   D (electrodes x time [20 samples before and 50 samples after start of phoneme] x phoneme trials) -- multiple trials of the same phoneme can be present, so this is just all trials concatenated together
%      and the variable [L] will identify the phoneme for each trial.)
%   L (phonemes trials x 1) -- this is the number labeling the phoneme in the
%                              matrix D
%   phnsub (cell array with unique(L) entries) -- these are the actual
%                                                 phoneme labels for L
% 
% For example:
% D = 256 electrodes x 70 time points x 7000 trials
% L = [1 1 1 ..... 33 33 33] (length 7000)
% phnsub = {'aa','ao',...'z'} (length 33, or however many phonemes).
% 

ecog = Phn.ecog_z;
channels = [17 25 33 41]; % Channels of interest
ref_chans = [18 26 34 42];
%ecog = ecog(channels,:,:) - ecog(ref_chans,:,:);

%% Remove hash #
is_hash = strcmpi(Phn.phn_names,'h#');
Phn.phn_names(is_hash) = [];
ecog(:,:,is_hash)= [];
onset = round(0.5*fs)-20;
offset = onset+169;
D = ecog(:,onset:offset,:);
%% Make phoneme name scheme;
names = unique(Phn.phn_names);
name_sizes = zeros(size(names));
for i = 1:length(names)
    name_sizes(i) = sum(strcmpi(Phn.phn_names,names(i)));
end
[names_sort, order] = sort(name_sizes,'descend');
names = names(order(1:33)); % Magic number?

names =  {  'aa'    'ao'    'ow'    'ah'    'ax' ...
    'uh'    'ux'    'uw'    'iy'    'ih' ...
    'ix'    'ey'    'eh'    'ae'    'aw' ...
    'ay'    'w'    'y'    'l'    'r' ...
    'dh'    'th'    'f'    's'    'sh' ...
     'z'    'v'    'p'    't'    'k'    'b' ...
     'd'    'g'    'm'    'n'    'ng'};
L = zeros(length(Phn.phn_names),1);
for i = 1:length(names)
    L = L + i.*double(strcmpi(Phn.phn_names, names(i)));
end
is_excluded = (L==0);
L(is_excluded) = [];
D(:,:,is_excluded) =[];

% size_lim = 2500;
% L = L(1:size_lim);
% D = D(:,:,1:size_lim);

phnsub = names;


% Create a window centered at 150 ms after phoneme onset (these bins are 10
% ms wide) -- this is hard coded and assumes 20 samples before phoneme
% onset is how data are stored, so you will need to change these
% appropriately depending on how you store your phoneme matrix
bef = 20;
win = 0.05*fs; % 50 ms window
win_center = 0.15;
tstart = floor(win_center*fs + bef - win/2);
tend = ceil(win_center*fs + bef + win/2);
tps = tstart:tend; % Time points to use in analysis

% run the stats
Lbs = phnsub; % Labels for 33 phonemes (renaming because this is what psicluster expects... kind of dumb)
Dp = squeeze(mean(D(:,tps,:),2)); % Take the mean phoneme response during this time period

tic;
% Preallocate matrices for ranksum tests
p_valp = zeros(length(phnsub), length(phnsub), size(Dp,1)); % matrix of p values
t_val = zeros(length(phnsub), length(phnsub), size(Dp,1)); % matrix of rank sum statistic values
h_valp = zeros(length(phnsub), length(phnsub), size(Dp,1)); % matrix of h values
for phn1 = 1:length(phnsub) % for each phoneme
    %fprintf(1, 'phoneme %d: %s\n', phn1, phnsub{phn1});
    for phn2 = 1:length(phnsub) % compare to every other phoneme
        for chan = 1:size(Dp,1) % for each electrode
            [p, h, tmp2] = ranksum(Dp(chan,L==phn2), Dp(chan,L==phn1), 'alpha', 0.01);
            p_valp(phn2,phn1,chan) = p;
            h_valp(phn2,phn1,chan) = h;
            t_val(phn2,phn1,chan) = tmp2.ranksum;
            if median(Dp(chan,L==phn1))>median(Dp(chan,L==phn2))
                h_valp(phn2,phn1,chan)=0; % make it a one-tailed test
            end
        end
    end
end

toc

smap = squeeze(sum(h_valp,2)); % This is the PSI map
figure; imagesc(smap)
ax = gca;
ax.XTick = 1:size(smap,2);
for i = 1:length(channels)
    ax.XTickLabel{i} = ['Ch ', num2str(channels(i)), ' - Ch ', num2str(ref_chans(i))];
end
%xlabel('Channel')
ax.YTick = 1:length(names);
ax.YTickLabel = names;
title('PSI for Bipolar Grid Electrodes');



ind = 1:size(D,1); % These are the electrode numbers

% Save to file
% outfile = sprintf('%s/%s/%s_PSI_%dto%d.mat', dpath, subject, subject, tstart,tend);
% fprintf(1,'Saving file %s\n', outfile);
% save(outfile,'smap','h_valp','p_valp','t_val','Lbs','ind');
% fprintf(1,'Done.\n');

% Can run this line to cluster the PSIs, or may want to do this after
% concatenating smap across subjects
%psi_thresh = 2;
%[ord_elec,ord_phone,tstsmap, t_elec, t_phone, z_elec, z_phone, good_chans] = psicluster(smap, Lbs, psithresh);
