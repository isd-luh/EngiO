% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function fig=plotParetoDesign(x_pareto, f_pareto, mSamples)
% PLOTPARETODESIGN creates plot of sampling points in design variable
% space for 2- and 3-dimensional problems. Non-dominated solutions are
% highlighted.
%
% INPUT:
% - x_pareto (m by n matrix of doubles)
%   Pareto-optimal design variables, where m is number of non-dominated
%   solutions and n number of design variables per solution
% - f_pareto (m by n matrix of doubles)
%   Pareto-optimal design variables, where m is number of non-dominated
%   solutions and n number of design variables per solution
% - mSamples (m by n matrix of doubles)
%   Sampled design varibles, where m is number of samples
%   and n number of design variables per sample
%
% OUPUT:
% - fig (figure object)
%   figure object with plot properties

fig=figure('name','Design Variable Space');
% Two dimensions
if size(mSamples, 2) == 2
    scatter(mSamples(:, 1), mSamples(:, 2), '.k')
    hold on
    scatter(x_pareto(:, 1), x_pareto(:, 2), 3, f_pareto(:, 1))
    
    if exist ('OCTAVE_VERSION', 'builtin')
        colorbar;
    else
        c = colorbar;
        c.Label.String = 'Objective Value 1';
    end
        
    colormap hsv
    
%     title ('Design Variable Space')
    xlabel('Design Variable 1')
    ylabel('Design Variable 2')
    legend('Sampled Points', 'Non-dominated Points')
end

% Three dimensions
if size(mSamples, 2) == 3
    extraSamples = setdiff(mSamples, x_pareto, 'rows');
    scatter3(extraSamples(:, 1), extraSamples(:, 2), extraSamples(:, 3), '.k', 'MarkerEdgeAlpha', 0.2)
    hold on
    scatter3(x_pareto(:, 1), x_pareto(:, 2), x_pareto(:, 3), 3, f_pareto(:, 1))
    
    if exist ('OCTAVE_VERSION', 'builtin')
        colorbar;
    else
        c = colorbar;
        c.Label.String = 'Objective Value 1';
    end
    
    colormap hsv
    
%     title ('Design Variable Space')
    xlabel('Design Variable 1')
    ylabel('Design Variable 2')
    zlabel('Design Variable 3')
    legend('Sampled points', 'Non-dominated points', 'Location', 'northeast')
end

end