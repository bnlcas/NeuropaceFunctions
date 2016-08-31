function time_range_bad = contains_badTimes(bad_Segments, time_range); % employs helper function - contains_badTimes
%% Takes a nx2 matrix where each row is an interval of bad time, and 
% a 1x2 matrix demaracting a region of interest in returns true iff the ROI
% overlaps any of the 1x2 rois in the badTimeSegments matrix
time_range_bad = false;
for i = 1:size(bad_Segments,1)
    bad_segment = bad_Segments(i,:);  
% % %     w1 = abs(diff(bad_segment)); % span of the bad time segment
% % %     w2 = abs(diff(time_range)); % span of the roi
% % %     if (max(bad_segment(2),time_range(2)) - min(bad_segment(1),time_range(1))) < (w1+w2);
% % %         time_range_bad = true
% % %     end

    if min(bad_segment(2),time_range(2)) >= max(bad_segment(1), time_range(1))
        time_range_bad = true;
    end

end

end
