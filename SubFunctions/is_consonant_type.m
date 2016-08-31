function [is_type] = is_consonant_type(Phn, list, varargin);
%% Create a boolean for whether phomemes match a list
%% as an optional input, speficy whether to use full comparison to the list or only look for a single character
is_type = false(size(Phn.phn_names));
for i = 1:length(list)
    test = list{i};
    if ~isempty(varargin)
        if varargin{1}
            is_type = is_type | substrcmp(Phn.phn_names,list{i});
        else
            is_type = is_type | strcmpi(Phn.phn_names,list{i});
        end
    else
        is_type = is_type | strcmpi(Phn.phn_names,list{i});
    end
end

