function phn_table = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   PHN_TABLE = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   PHN_TABLE = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   phn_table = importfile('fadg0_si1279.phn', 1, 24);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/11/06 15:02:17

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
dataArray(4) = [];

fclose(fileID);

%% Output
phn_table = table(dataArray{1:3}, 'VariableNames', {'start','stop','phn'});



