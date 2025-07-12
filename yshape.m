function plotYShape(theta_deg)
% PLOTYSHAPE Plot a single Y-shape with specified angle
%   plotYShape(theta_deg) plots a Y-shape with given angle in degrees
%   - Segments are colored red
%   - Points are colored green

    % Default parameters
    L = 125; % mm (length of main stem)
    
    % Convert angle to radians
    theta = deg2rad(30);
    
    % Calculate coordinates
    P1 = [0, 0];                    % Base point
    P2 = [0, L];                    % Vertical top
    P3 = [L*cos(theta), L+L*sin(theta)];  % Right branch
    P4 = [-L*cos(theta), L+L*sin(theta)]; % Left branch
    
    % Create figure
    figure;
    hold on;
    axis equal;
    grid on;
    
    % Plot the Y-shape segments in red
    plot([P1(1), P2(1)], [P1(2), P2(2)], 'r-', 'LineWidth', 2.5); % Vertical stem
    plot([P2(1), P3(1)], [P2(2), P3(2)], 'r-', 'LineWidth', 2.5); % Right branch
    plot([P2(1), P4(1)], [P2(2), P4(2)], 'r-', 'LineWidth', 2.5); % Left branch
    
    % Mark vertices in green
    scatter([P1(1), P2(1), P3(1), P4(1)], [P1(2), P2(2), P3(2), P4(2)], ...
            100, 'g', 'filled', 'MarkerEdgeColor', 'k');
    
    % Label points
    text(P1(1), P1(2)-15, 'P1 (0,0)', 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(P2(1), P2(2)+15, sprintf('P2 (0,%d)', L), 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(P3(1)+10, P3(2), sprintf('P3 (%.1f,%.1f)', P3(1), P3(2)), 'FontSize', 10);
    text(P4(1)-10, P4(2), sprintf('P4 (%.1f,%.1f)', P4(1), P4(2)), ...
        'HorizontalAlignment', 'right', 'FontSize', 10);
    
    % Format plot
    title(sprintf('Y-Shape (θ = %d°)', theta_deg), 'FontSize', 12);
    xlabel('X coordinate (mm)', 'FontSize', 11);
    ylabel('Y coordinate (mm)', 'FontSize', 11);
    
    % Set axis limits
    axis_buffer = 50;
    xlim([P4(1)-axis_buffer, P3(1)+axis_buffer]);
    ylim([-axis_buffer, P3(2)+axis_buffer]);
    
    hold off;
end