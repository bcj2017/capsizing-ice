function drawCurvedArrow(theta_start, theta_end, tri_size, center_x, center_y, theta_arrowhead, arrowcolor)
    % Define the circle segment
    theta = linspace(theta_start, theta_end, 100); % Generate points for smooth curve
    base_center_x = center_x;
    base_center_y = center_y;
    
    % Calculate the height of the equilateral triangle
    height = sqrt(3) / 2 * tri_size;
    
    % Define the vertices of the equilateral triangle
    v1 = [0, height]; % The top vertex of the triangle
    v2 = [-tri_size / 2, 0]; % Bottom left vertex
    v3 = [tri_size / 2, 0]; % Bottom right vertex

    % Rotate the triangle vertices
    R = [cos(theta_arrowhead), -sin(theta_arrowhead); sin(theta_arrowhead), cos(theta_arrowhead)];
    v1_rot = (R * v1')' + [base_center_x, base_center_y];
    v2_rot = (R * v2')' + [base_center_x, base_center_y];
    v3_rot = (R * v3')' + [base_center_x, base_center_y];
    
    % Plot the filled equilateral triangle
    fill([v1_rot(1), v2_rot(1), v3_rot(1)], [v1_rot(2), v2_rot(2), v3_rot(2)], arrowcolor, 'EdgeColor', 'none', 'LineWidth', 2);
end
