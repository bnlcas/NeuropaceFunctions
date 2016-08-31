function [evnt_redir] = transform_evnt_dir(evnt, root_dir)
%% takes a TIMIT evnt structured array and converts the data path to be realligned to the root directory
% of the local disc

evnt_redir = evnt;

for i = 1:length(evnt)
    path_original = strsplit(evnt(i).dpath,'/');
    if strcmpi(path_original(end-2), 'EC71')
        path_tail = path_original((end-3):end);
    else
        path_tail = path_original((end-2):end); % the last bit of the path (/data/Session/Block)
    end
    
    evnt_redir(i).dpath = [root_dir, '/', strjoin(path_tail,'/')];
    
    timit_dir_original = strsplit(evnt(i).wname,'/');
    path_tail = timit_dir_original((end-3):end);
    evnt_redir(i).wname = [root_dir, '/' strjoin(path_tail,'/')];
end

end