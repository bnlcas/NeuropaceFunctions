function [is_bad_trial] = find_bad_trials_evnt(evnt)
%% This function reads through the evnt structed array for a block of RNS recording
% and find trials that occur during bad time segments

is_bad_trial = false(length(evnt),1);

%% Handle EC71 (no offset on raw / timelocked ecog)
if substrcmp(evnt(1).expt,'EC71');
    for i = 1:length(evnt);
        dpath = evnt(i).dpath;
        if i == 248
            a = 1;
        end
        if isempty(evnt(i).StartTime)
            is_bad_trial(i) = true;
        else
            load([dpath '/Artifacts/badTimeSegments.mat'])
        
            evnt_time_range = [evnt(i).StartTime evnt(i).StopTime];
            is_bad_trial(i) = contains_badTimes(badTimeSegments, evnt_time_range);
        end
    end
else % RNS DATA:
    for i = 1:length(evnt)
        dpath = evnt(i).dpath;
        tmp = strsplit(dpath,'/');
        tmp{end-1} = [tmp{end-1} '_RAW'];
        dpath_artifacts = strjoin(tmp,'/');
        load([dpath_artifacts '/Artifacts/badTimeSegments.mat'])

        evnt_time_range = [evnt(i).StartTime evnt(i).StopTime];
        is_bad_trial(i) = contains_badTimes(badTimeSegments, evnt_time_range);
    end   
end 


