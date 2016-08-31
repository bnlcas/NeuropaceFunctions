function ecog_prcd = remove_64_samp_artifact(ecog_ch);
%% Finds the starting index of a periodic distortion 
% occuring every 64 samples and smooths it out


%% Find index
%delta = diff(ecog_ch);
glitch_rate = 64;
test = zeros(length(glitch_rate),1);

% for i = 1:(glitch_rate-1)
%     test(i) = mean(delta(i:glitch_rate:length(delta)));
%     error(i) = std(delta(i:glitch_rate:length(delta)));
% end
% [shift,glitch_ind] = max(abs(test));

for i = (1:glitch_rate)
     tmp(i) = mean(ecog_ch(i:glitch_rate:length(ecog_ch)));
end
% glitch results in values being distinct from the average
[shift, glitch_ind] = max(abs(tmp-mean(tmp)));

%% Smooth over glitched indecies:
ecog_prcd = ecog_ch;
%test = ecog_ch;
glitched_inds = glitch_ind:glitch_rate:length(ecog_ch);
if glitched_inds(1) < 3 
    glitched_inds = (glitch_ind+glitch_rate):glitch_rate:length(ecog_ch);
end
if glitched_inds(end) > (length(ecog_ch)-2)
    glitched_inds = (glitch_ind):glitch_rate:(length(ecog_ch)-glitch_rate);
end

for i = glitched_inds
    x_inds = [-2:-1, 1:2];
    neighbors = ecog_ch(i+x_inds);
    y_spline = spline(x_inds, neighbors, -2:2);
    ecog_prcd(i) = y_spline(3);
    %test(i) = ecog_ch(i) + shift;
end


end
