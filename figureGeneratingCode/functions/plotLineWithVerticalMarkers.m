function plotLineWithVerticalMarkers(pt1, pt2, height, linewidth)
    % Plot the main line
    hold on;
    plot([pt1(1), pt2(1)], [pt1(2), pt2(2)], 'LineWidth', linewidth, 'Color', 'k');
    
    % Calculate vertical lines' endpoints
    % Midpoint of the line segment
    midpoint1 = pt1;
    midpoint2 = pt2;
    
    % Vertical lines
    plot([midpoint1(1), midpoint1(1)], [midpoint1(2) - height/2, midpoint1(2) + height/2], ...
         'LineWidth', linewidth, 'Color', 'k');
    plot([midpoint2(1), midpoint2(1)], [midpoint2(2) - height/2, midpoint2(2) + height/2], ...
         'LineWidth', linewidth, 'Color', 'k');
end
