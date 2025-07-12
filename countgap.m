clc;
clear;
close all;

% Parameters
omeg_mult = 10;
angles = 5:70;  % Angles to process
output_file = 'bandgaps_summary_with_width.csv';

% Initialize output data
all_bandgaps = [];

% Process each angle
for angle = angles
    try
        % Load data for current angle
        eigenfrequency_data = load(sprintf('Lin_Norm_freq30.dat', angle));
        wave_number_data = load('Wave_num.dat')';
        
        % Identify band gaps
        sdat = unique(sort(eigenfrequency_data(:)));  % Flatten, sort, and get unique frequencies
        gap_indices = find(diff(sdat) >= 1e-1);       % Find indices where a significant gap exists
        
        % Store bandgap information
        if ~isempty(gap_indices)
            for i = 1:length(gap_indices)
                bandgap_start = sdat(gap_indices(i));
                bandgap_end = sdat(gap_indices(i) + 1);
                bandgap_width = bandgap_end - bandgap_start;
                
                % Add to output array
                all_bandgaps = [all_bandgaps; angle, i, bandgap_start, bandgap_end, bandgap_width];
            end
        else
            % If no bandgap found, record angle with NaN values
            all_bandgaps = [all_bandgaps; angle, NaN, NaN, NaN, NaN];
        end
        
    catch ME
        warning('Failed to process angle %d: %s', angle, ME.message);
        all_bandgaps = [all_bandgaps; angle, NaN, NaN, NaN, NaN];
    end
end

% Write results to CSV file
header = {'Angle', 'Bandgap_Number', 'Start_Frequency', 'End_Frequency', 'Bandgap_Width'};
fid = fopen(output_file, 'w');
fprintf(fid, '%s,%s,%s,%s,%s\n', header{:});
fclose(fid);
dlmwrite(output_file, all_bandgaps, '-append', 'precision', '%.6f', 'delimiter', ',');

disp(['Bandgap data with width calculations saved to ' output_file]);

% Optional: Display summary statistics
if ~all(isnan(all_bandgaps(:,5)))
    fprintf('\nBandgap width statistics:\n');
    fprintf('Maximum bandgap width: %.4f\n', max(all_bandgaps(:,5)));
    fprintf('Minimum bandgap width: %.4f\n', min(all_bandgaps(:,5)));
    fprintf('Average bandgap width: %.4f\n', mean(all_bandgaps(:,5), 'omitnan'));
end