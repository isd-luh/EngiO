% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

% EXAMPLECON demonstrates usage for constrained optimization
% problems

% Add folders to path
addpath('OptAlgorithms', 'TestFunctionsCon', 'Penalties',...
    'Plots', '-frozen');

% Set optimization problem
prob = @Simionescu;

% Set lower and upper bounds for problem
x_lb = [-1.25, -1.25];
x_ub = [1.25 ,  1.25];

% Create instance of optimization algorithm
algo = GlobalPattern();

% Create constraint handler
pen = PenalizedObjectiveFunction(prob, {ExteriorLinearPenalty(.1)}, {});

% Settings for optimization run
aOptions = struct('maxEvals', 1000, 'saveStates', false);
aParams  = struct('nTrack' ,   20);

% Run optimizer
[x_opt, y_opt, nEvals, vSamples, vResults] = algo.optimize(@(x) pen.evaluate(x), [], x_lb, x_ub, aOptions, aParams);

% Print result
result = prob(x_opt);
disp(['x = [', num2str(x_opt(1)), ', ', num2str(x_opt(2)), ']; f = ', num2str(result{1}), '; g = ', num2str(result{2})])

% Create plots
fig(1)=plotObjectiveHistory(vSamples, vResults, prob);
fig(2)=plotConstraintHistory(vSamples, vResults, prob);

% Remove folders from path
rmpath('OptAlgorithms', 'TestFunctionsCon', 'Penalties', 'Plots');