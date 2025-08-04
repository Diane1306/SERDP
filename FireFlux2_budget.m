clc
clear
%%
% Main script
data_dir = '/Users/diane_wt/Downloads/work/FireFlux2/';

% List of text files to load
file_names = {'m6m_10hz.txt', 'm10m_10hz.txt', 'm20m_10hz.txt'};

% Loop through each file and load the data
for i = 1:length(file_names)
    % Construct the full file path
    file_path = fullfile(data_dir, file_names{i});

    % Load the data from the text file into a table
    % Assuming the first row contains variable names
    data = readtable(file_path);

    % Extract height from the file name (e.g., 'm06m' -> '6m')
    height = extractBetween(file_names{i}, 'm', '_10hz');

    % Create new variable names based on the height
    new_var_names = strcat(data.Properties.VariableNames, '_', height);

    % Rename the variables in the table
    data.Properties.VariableNames = new_var_names;

    % Assign the modified table to a variable in the workspace
    assignin('base', strcat('data_', height{1}), data);
end