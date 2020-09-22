% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

% BENCHMARKSINGLEOBJ benchmarks optimzers with single-objective
% unconstrained test functions

% Octave: requires the statistics package

% Add folders to path
addpath('OptAlgorithms', 'TestFunctionsUncon', '-frozen')

% Set optimization problems
aProb = {
    struct('name', 'Ackley', 'fun', @Ackley, 'vartype', [], 'x_lb', [-5, -5], 'x_ub', [5, 5]),...
    struct('name', 'Beale', 'fun', @Beale, 'vartype', [],  'x_lb', [-4.5 -4.5], 'x_ub', [4.5 4.5]),...
    struct('name', 'Camel 6', 'fun', @Camel6, 'vartype', [],  'x_lb', [-5 -5], 'x_ub', [5 5]),...
    struct('name', 'Cross in Tray', 'fun', @Crossintray, 'vartype', [],  'x_lb', [-10, -10], 'x_ub', [10, 10]),...
    struct('name', 'Eggholder', 'fun', @Eggholder, 'vartype', [],  'x_lb', [-512, -512], 'x_ub', [512, 512]),...
    struct('name', 'Himmelblau', 'fun', @Himmelblau, 'vartype', [],  'x_lb', [-6, -6], 'x_ub', [6, 6]),...
    struct('name', 'Rastrigin', 'fun', @Rastrigin, 'vartype', [],  'x_lb', [-5.12, -5.12], 'x_ub', [5.12, 5.12]),...
    struct('name', 'Rosenbrock', 'fun', @Rosenbrock, 'vartype', [],  'x_lb', [-2, -1], 'x_ub', [2, 3]),...
    struct('name', 'Schwefel', 'fun', @Schwefel, 'vartype', [],  'x_lb', [-500, -500], 'x_ub', [500, 500]),...
    struct('name', 'Sphere', 'fun', @Sphere, 'vartype', [],  'x_lb', [-5.12, -5.12], 'x_ub', [5.12, 5.12]),...
    };

% Settings for optimization run
aOptions = struct('outputStatus', false);

plotstyle = 'contour'; % or 'spatial';

% Loop over optimization problems
for iProb=1:size(aProb, 2)
    prob = aProb{iProb};
    
    fprintf('\n-- Problem: %s --\n', prob.name)
    fprintf('Optimizer\t\tEvaluations\t\tBest Objective\n')
    
    % Create instances of optimization algorithm
    vAlgo = {
        Random(), Grid(), GlobalPattern(),...
        ParticleSwarm(), EvolutionStrategy(), Harmony(),...
        Genetic(), SimulatedAnnealing()};
    
    % Prepare contour plot
    [X, Y] = meshgrid(prob.x_lb(1):.01*(prob.x_ub(1) - prob.x_lb(1)):prob.x_ub(1),...
        prob.x_lb(2):.01*(prob.x_ub(2) - prob.x_lb(2)):prob.x_ub(2));
    Z = zeros(size(X, 1), size(X, 2));
    logZ = Z;
    
    if length(aProb{iProb}.x_lb)==2
        for i = 1:size(X, 1)
            for j = 1:size(X, 2)
                x = [X(1, i), Y(j, 1)];
                Z(j, i) = prob.fun(x);
            end
        end
        logZ = log(Z-min(min(Z)));
        
        figure('name',prob.name);
    end
    
    numAlgo=0;
    
    % Loop over optimization algorithms
    for algo = vAlgo
        strAlgo = class(algo{1});
        
        % Run optimization
        [x_opt, y_opt, nEvals, mSamples, vResults] = algo{1}.optimize(prob.fun, prob.vartype, prob.x_lb, prob.x_ub, aOptions);
        
        % Evaluate objective function
        fObjective = prob.fun(x_opt(1,:));
        fprintf('%-20s%d\t\t%e\n', strAlgo, nEvals, fObjective)
        
        % mSamples(n, 3) >> col. 1 - 2 with design values | col. 3 with function values
        if length(aProb{iProb}.x_lb)==2
            numAlgo=numAlgo+1;
            
            subplot(ceil(numel(vAlgo)/4), 4, numAlgo)
            strTitle = sprintf(strAlgo);
            
            % Extract best sample coordinates
            [vfSortedMins, viSortedIndices] = sort(vResults);
            
            % Look for clusters in the sample data
            sampPercentile = .05; % percentile of samples taken into account
            cluThresPercentile = .3; % percentile of euclidean distances in clusters
            maxClu = 6; % maximum number of cluster
            minSamp = 10; % minimum number of samples for cluster analysis to be useful
            
            % Filter samples by percentile
            ind = viSortedIndices(1 : round(sampPercentile * numel(vResults)) );
            
            % Create set of samples near the optimum residual value
            chosen = mSamples(ind, :);
            
            % Get number of clusters by hierarchical clustering
            cluDist = pdist(chosen, 'euclidean');
            cluZ = linkage(cluDist);
            maxCluInconsistency = max(cluZ(:, 3) );
            cluThres = maxCluInconsistency * cluThresPercentile;
            mask = cluZ(:, 3) > cluThres;
            nClu = sum(mask) + 1;
            
            % If number of clusters is too large, there usually is just one
            if nClu > maxClu
                nClu = 1;
            end
            
            % Get cluster centers using kmeans
            [idx, C] = kmeans(chosen, nClu);
            vfOpt = mSamples(viSortedIndices(1), :);
            
            if strcmp(plotstyle, 'contour')
                contour(X, Y, logZ, 13,'HandleVisibility','off');
                title(strTitle);
                hold on
                plot(mSamples(:, 1), mSamples(:, 2), '.k', 'MarkerSize', 8);
                plot(vfOpt(1), vfOpt(2), '.r', 'MarkerSize', 20);
                
                plot(C(:, 1), C(:, 2), '.g', 'MarkerSize', 10);
               
                axis equal
                xlabel('Design Variable 1')
                ylabel('Design Variable 2')

            elseif strcmp(plotstyle, 'spatial')
                surf(X, Y, Z);
                title(strTitle);
                hold on
                
                plot3(vfOpt(1), vfOpt(2), vResults(viSortedIndices(1)), '.y', 'MarkerSize', 20);
                plot3(mSamples(:, 1), mSamples(:, 2), vResults, '.m', 'MarkerSize', 8);
            end
            hold off 
        end
    end
     legend({'Samples','Best Sample','Minimum'});
end


% Remove folders from path
rmpath('OptAlgorithms', 'TestFunctionsUncon')