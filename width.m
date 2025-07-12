clc;
clear;
close all;

% Parameters
frequency_range = 5:70;  % Frequencies from 5 to 70 Hz
output_filename = 'bandgap_with_only_width.csv';
gap_threshold = 1e-1;   % Minimum gap size to be considered a bandgap

% Initialize dataset with only required columns
headers = {'Angle', 'Bandgap_Number', 'Width'};
dataset = cell(0, length(headers));

% Process each frequency file
for freq = frequency_range
    % Load data
    filename = sprintf('Lin_Norm_freq%d.dat', freq);
    
    if ~exist(filename, 'file')
        continue;  % Skip silently
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
            bandgap_width = bandgap_end - bandgap_start;
            
            % Use frequency as Angle (modify if you have specific angle values)
            angle = freq;  % Replace with actual angle if available
            
            % Add to dataset with only required columns
            dataset(end+1, :) = {angle, gap_num, bandgap_width};
        end
    end
end

% Save as CSV
if ~isempty(dataset)
    tbl = cell2table(dataset, 'VariableNames', headers);
    writetable(tbl, output_filename);
end