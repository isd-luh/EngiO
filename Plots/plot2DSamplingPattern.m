% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function fig=plot2DSamplingPattern(x_opt, mSamples, prob, x_lb, x_ub)
% PLOT2DSAMPLINGPATTERN creates a contour plot of the objective function and
% marks sampling points and best solution.
%
% INPUT:
% - x_opt (vector of doubles)
%   Best design variables
% - mSamples (m by n matrix of doubles)
%   Sampled design varibles, where m is number of samples
%   and n number of design variables per sample
% - prob (function handle)
%   Objective function
% - x_lb (vector of doubles)
%   Lower boundaries of design variables
% - x_ub (vector of doubles)
%   Upper boundaries of design variables
%
% OUPUT:
% - fig (figure object)
%   figure object with plot properties

% Prepare objective function for contour plot
[X, Y] = meshgrid(x_lb(1):.01*(x_ub(1) - x_lb(1)):x_ub(1), x_lb(2):.01*(x_ub(2) - x_lb(2)):x_ub(2));
Z = zeros(size(X, 1), size(X, 2));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x = [X(1, i), Y(j, 1)];
        Z(j, i) = prob(x);
    end
end

% Plot objective function
fig=figure('name','Sampling Pattern');
if isempty(find(Z<=0))
    contour(X,Y,log(Z), 15,'HandleVisibility','off');
else
    contour(X,Y,Z, 15,'HandleVisibility','off');
end
hold on
% plot(mSamples(:, 1), mSamples(:, 2), '.k', 'MarkerSize', 8);
scatter(mSamples(:,1),mSamples(:,2), '.k');
plot(x_opt(1), x_opt(2), 'or','MarkerSize', 5);
grid on
axis equal
xlabel('Design Variable 1')
ylabel('Design Variable 2')
legend({'Samples','Best Sample'});

end