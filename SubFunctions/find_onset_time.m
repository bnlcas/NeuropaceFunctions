function start_times = find_onset_time(Audio, fs, latency)
%% Takes an audio input and finds time points when sound pulses start

start_times = zeros(size(Audio));

power = Audio.^2;

thresh = 0.5*std(power);
exceedes_thresh = (power > thresh);

i = 1;
latency = fs*latency;     % corresponds to a 1 second latency period

while i < length(Audio)
    if exceedes_thresh(i)
        start_times(i) = 1;
        i = i + latency;
    end
    i = i+1;
end

end
    
