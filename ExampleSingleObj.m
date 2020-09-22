% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

% EXAMPLESINGLEOBJ demonstrates usage for single-objective optimization
% problems

% Add folders to path
addpath('OptAlgorithms', 'TestFunctionsUncon', 'Plots', '-frozen');

% Set optimization problem
prob = @Schwefel;

% Set lower and upper bounds for problem
x_lb = [-500, -500];
x_ub = [ 500,  500];

% Create instance of optimization algorithm
algo = Genetic();

% Settings for optimization run
aOptions = struct('maxEvals', 1000, 'saveStates', false);
aParams  = struct('popSize' ,   30);

% Run optimizer
[x_opt, y_opt, nEvals, mSamples, vResults]...
    = algo.optimize(prob, [], x_lb, x_ub, aOptions, aParams);

vResults = vResults - min(vResults);

% Retrieve objective value for optimum
y = prob(x_opt);
fprintf('\nOptimal Design Variables:\t(%f, %f)\nOptimal Objective Value:\t%e\n', x_opt(1), x_opt(2), y);

% Create plots
fig(1)=plot2DSamplingPattern(x_opt, mSamples, prob, x_lb, x_ub);
fig(2)=plotDesignHistory(mSamples);
fig(3)=plotObjectiveHistory(mSamples, vResults, prob);

% Remove folders from path
rmpath('OptAlgorithms', 'TestFunctionsUncon', 'Plots');