function output = substrcmp(input, string)
%% Takes an input 'input' of string or cell array of strings and returns a boolean if it contains the input 'string'

%%Convert a single string input to a cell array of string with length 1
if isstr(input)
    input = {input};
end
num_elements = length(input);
output = false(num_elements,1);

%% Create array of index at which string occurs in each element of input
contains_string = strfind(input, string);

%% For nonemplty elements, return 0
for i = 1:num_elements
    element = cell2mat(contains_string(i));
    if ~isempty(element)
        output(i) = true;
    end
end

end
    