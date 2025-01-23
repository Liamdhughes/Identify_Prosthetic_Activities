% Copyright, M.Bencsik, M.Bisele L.D.Hughes, 2025

function y = graph_scores_5_polygon(scores, DFA_scores, A, B, boundaries)
    % Calculate overall min and max for x and y
    overall_min_x = min(DFA_scores(:,1));
    overall_max_x = max(DFA_scores(:,1));
    overall_min_y = min(DFA_scores(:,2));
    overall_max_y = max(DFA_scores(:,2));
    
    % Calculate center points
    center_x = (overall_max_x + overall_min_x) / 2;
    center_y = (overall_max_y + overall_min_y) / 2;
    
    % Define axis limits to ensure each plot has a range of 10 units
    xlims = [center_x - 30, center_x + 30];
    ylims = [center_y - 30, center_y + 30];
    
    % Initialize plot parameters
    param = zeros(5,2);
    param(:,1) = [1; boundaries(1)+1; boundaries(2)+1; boundaries(3)+1; boundaries(4)+1];
    param(:,2) = [boundaries(1); boundaries(2); boundaries(3); boundaries(4); boundaries(5)];

    % Plot and hold to add more plots
    hold on;
    colors = {'r', 'k', 'c', 'g', 'b'}; % Colors for each polygon
    
    % Loop through each boundary set and plot
    for i = 1:size(param, 1)
        DFA_subset = DFA_scores(param(i,1):param(i,2), :);
        plot(DFA_subset(:,1), DFA_subset(:,2), 'o', 'MarkerEdgeColor', colors{i}, 'MarkerFaceColor', colors{i});
        
        % Calculate and plot the boundary for each subset
        polygon_boundary = boundary(DFA_subset(:,1), DFA_subset(:,2), 0.1);
        plot(DFA_subset(polygon_boundary,1), DFA_subset(polygon_boundary,2), 'Color', colors{i}, 'LineWidth', 1.5);
    end

    % Set common axis limits after all plots
    axis([xlims ylims]);

    % Additional plot formatting
    grid on;
    xlabel('DF score No 1');
    ylabel('DF score No 2');
    legend('Level', '', 'Incline', '', 'Decline', '', 'Upstairs', '', 'Downstairs', '');
    legend('Location', 'northeastoutside');
    hold off; % Release the hold to allow for new plots
end
