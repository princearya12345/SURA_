clc;
clear;
close all;

% Parameters
omeg_mult = 10;

% Load data
eigenfrequency_data = load('Lin_Norm_freq30.dat'); % Eigenfrequencies (y-axis)
wave_number_data = load('Wave_num.dat')';      % Wave numbers (x-axis)

% Plot the dispersion diagram
figure;
hold on;
plot(wave_number_data, eigenfrequency_data, '-b', 'LineWidth', 2);
ylim([0, omeg_mult]);

% Identify band gaps
sdat = unique(eigenfrequency_data(:));        % Flatten and get unique frequencies
gap_indices = find(diff(sdat) >= 1e-1);        % Find indices where a significant gap exists

% Highlight all band gaps
if ~isempty(gap_indices)
    disp('Bandgaps found:');
    for i = 1:length(gap_indices)
        % Calculate start and end of each band gap
        bandgap_start = sdat(gap_indices(i));
        bandgap_end = sdat(gap_indices(i) + 1);
        
        % Draw each band gap as a horizontal rectangle
        rectangle('Position', [min(wave_number_data), bandgap_start, ...
                  max(wave_number_data) - min(wave_number_data), bandgap_end - bandgap_start], ...
                  'FaceColor', 'k', 'EdgeColor', 'none');%, 'FaceAlpha', 0.3); % Semi-transparent
        
        % Display bandgap details in the command window
        disp(['Bandgap ', num2str(i), ': [', num2str(bandgap_start), ', ', num2str(bandgap_end), ']']);
    end
else
    disp('No bandgaps found.');
end

% Annotate plot
xlabel('Bloch wave vector κ');
ylabel('Normalized Frequency Ω');
set(gca, 'FontSize', 12); % Adjust font size for readability
hold off;