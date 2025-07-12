clc;
clear;
close all;

% Parameters
frequency_range = 5:70;  % Frequencies from 5 to 70 Hz
output_filename = 'dataset_with_start_end.csv';
gap_threshold = 1e-1;   % Minimum gap size to be considered a bandgap

% Initialize dataset
headers = {'Angle', 'Bandgap_Number', 'Start_Frequency', 'End_Frequency'};
dataset = cell(0, length(headers));

% Process each frequency file
for freq = frequency_range
    % Load data
    filename = sprintf('Lin_Norm_freq%d.dat', freq);
    
    if ~exist(filename, 'file')
        warning('File %s not found. Skipping.', filename);
        continue;
    end
    
    eigenfrequency_data = load(filename);
    
    % Identify band gaps
    sdat = unique(eigenfrequency_data(:));  % Sorted unique frequencies
    gap_indices = find(diff(sdat) >= gap_threshold);
    
    % Store all bandgaps for this frequency
    if ~isempty(gap_indices)
        for gap_num = 1:length(gap_indices)
            bandgap_start = sdat(gap_indices(gap_num));
            bandgap_end = sdat(gap_indices(gap_num) + 1);
            
            % Add to dataset
            dataset(end+1, :) = {freq, gap_num, bandgap_start, bandgap_end};
        end
    end
end

% Save as CSV
if ~isempty(dataset)
    tbl = cell2table(dataset, 'VariableNames', headers);
    writetable(tbl, output_filename);
    fprintf('Bandgap dataset saved to: %s\n', output_filename);
    
    % Display sample
    disp('Sample of the dataset:');
    disp(tbl(1:min(10,height(tbl)),:)); % Show first 10 rows or fewer
else
    disp('No bandgaps found in any frequency files.');
end