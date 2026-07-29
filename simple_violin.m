


function simple_violin(data, positions, color, width)
    % Simple violin plot using kernel density estimation
    if nargin < 4, width = 0.3; end
    if nargin < 3, color = [0.7 0.7 0.7]; end
    
    for i = 1:length(data)
        y = data{i};
        if isempty(y), continue; end
        
        % Kernel density estimation
        [f, xi] = ksdensity(y);
        f = f / max(f) * width; % Normalize width
        
        % Plot violin (mirrored density)
        patch([positions(i) - f, positions(i) + flip(f)], ...
              [xi, flip(xi)], color, 'EdgeColor', 'k', 'LineWidth', 1, 'FaceAlpha', 0.6);
        hold on;
        
        % Add median line
        plot([positions(i)-width/2, positions(i)+width/2], [median(y), median(y)], ...
             'k-', 'LineWidth', 2);

        % Add mean 
        % plot(positions(i), mean(y), 'rd', 'MarkerSize', 14, 'LineWidth', 2, 'MarkerFaceColor', 'r');
        
        % Add quartiles
        q = quantile(y, [0.25, 0.75]);
        plot([positions(i), positions(i)], q, 'k-', 'LineWidth', 3);
    end
end