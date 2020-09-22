% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function fig=plotParetoObjective(f_pareto, mResults)
% PLOTPARETODESIGN creates plot of sampling points in objective
% space for 2- and 3-objective problems. Non-dominated solutions are
% highlighted.
%
% INPUT:
% - f_pareto (m by n matrix of doubles)
%   Pareto-optimal design variables, where m is number of non-dominated
%   solutions and n number of objectives
% - mResults (m by n matrix of doubles)
%   Sampled objective values, where m is number of samples
%   and n number of objectives per sample
%
% OUPUT:
% - fig (figure object)
%   figure object with plot properties

fig=figure('name','Objective Value Space');

% Two objectives
if size(f_pareto, 2) == 2
    subplot(2,1,1);
    scatter(f_pareto(:, 1), f_pareto(:, 2),3,f_pareto(:, 1));
    hold on
    
    if exist ('OCTAVE_VERSION', 'builtin')
        colorbar;
    else
        c = colorbar;
        c.Label.String = 'Objective Value 1';
    end
    
    colormap hsv
    
    grid on
    
%     title ('Objective Value Space')
    xlabel('Objective Value 1')
    ylabel('Objective Value 2')
    legend('Non-dominated Points')
    
    subplot(2,1,2);
    scatter(mResults(:, 1), mResults(:, 2), '.k');
    hold on
    scatter(f_pareto(:, 1), f_pareto(:, 2), '.r');
    grid on
    
%     title ('Objective Value Space')
    xlabel('Objective Value 1')
    ylabel('Objective Value 2')
    legend('Sampled Points', 'Non-dominated Points')
end

% Three objectives
if size(f_pareto, 2) == 3
    subplot(2,1,1);
    scatter3(f_pareto(:, 1), f_pareto(:, 2), f_pareto(:, 3), '.k', 'MarkerEdgeAlpha', 0.2);
    hold on
    scatter3(ys(:, 1), ys(:, 2), ys(:, 3), 3, ys(:, 1));
    
    if exist ('OCTAVE_VERSION', 'builtin')
        colorbar;
    else
        c = colorbar;
        c.Label.String = 'Objective Value 1';
    end
    
    colormap hsv
    
    grid on
    
%     title ('Objective Value Space')
    xlabel('Objective Value 1')
    ylabel('Objective Value 2')
    zlabel('Objective Value 3')
    legend('Non-dominated Points', 'Points on Convex Front', 'Location', 'northeast')
    
    subplot(2,1,2);
    scatter3(mResults(:, 1), mResults(:, 2), mResults(:, 3), '.k');
    hold on
    scatter3(f_pareto(:, 1), f_pareto(:, 2), f_pareto(:, 3), '.r');
    grid on
    
%     title ('Objective Value Space')
    xlabel('Objective Value 1')
    ylabel('Objective Value 2')
    zlabel('Objective Value 3')
    legend('Sampled Points', 'Non-dominated Points', 'Location', 'northeast')
end

end