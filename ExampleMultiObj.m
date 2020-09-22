% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

% EXAMPLEMULTIOBJ demonstrates usage for multi-objective optimization
% problems

% Add folders to path
addpath('OptAlgorithms', 'TestFunctionsUncon', 'Plots', '-frozen');

% Set optimization problem
prob = @Poloni;

% Set lower and upper bounds for problem
x_lb = [-pi -pi];
x_ub = [ pi  pi];

% Get dimension of outputs
nObj = prob([]);

% Create instance of optimization algorithm
algo = GlobalPatternMO();

% Settings for optimization run
aOptions = struct('maxEvals', 10000, 'saveStates', false);
aParams = struct();

% Run optimizer
[x_pareto, f_pareto, nEvals, mSamples, mResults] = algo.optimize(prob, [], x_lb, x_ub, aOptions, aParams);

% Create plots
fig(1)=plotParetoDesign(x_pareto, f_pareto, mSamples);
fig(2)=plotParetoObjective(f_pareto, mResults);

% Remove folders from path
rmpath('OptAlgorithms', 'TestFunctionsUncon', 'Plots');

