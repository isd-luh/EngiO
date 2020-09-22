% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef GlobalPatternMO < Optimizer
    % GLOBALPATTERNMO is a multi-objective version of GLOBALPATTERN. It is
    % derived from Optimizer.
    %
    % LITERATURE:
    %   not yet published.
    
    methods(Access = public)
        function this = GlobalPatternMO()
            
            this.defaultparams = struct('nTrack',      10, ...                      % Local minima to track
                'nResolution',   2^11, ...                      % Resolution of the grid
                'nScalarizationResolution', 200,  ...           % Points sampled on the pareto front
                'vfNonlinearizations', 2.^linspace(0, 6, 7) ... % Exponents for nonlinear transforms
                );
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            nDims = numel(x_lb);
            
            this.states = struct('viResolution', zeros(nDims), ...  % Allocate space for resolution vector
                'miSampleCache', [], ...           % Cache for sample coordinates
                'vfObjectiveValues', [], ...       % Cache for objective values
                'nNonlinearizations', numel(this.params.vfNonlinearizations), ...
                'mfParetoMultipliers', [], ...     % Vectors for scalarization
                'viHallOfFame', [], ...            % Indices of samples currently on the front
                'viConvexFront', [] ...            % Samples on the convex front (elitists)
                );
            
            % Initialize resolution per dimension
            % Use given resolution for real dimensions
            this.states.viResolution = ones(1,nDims) * this.params.nResolution;
            % Use actual resolution for integer dimensions
            if ~isempty(vartype)
                for m = 1:nDims
                    if strcmpi(vartype{m}, 'real')
                        this.states.viResolution(m) = this.params.nResolution;
                    elseif strcmpi(vartype{m}, 'int')
                        this.states.viResolution(m) = x_ub(m) - x_lb(m);
                    end
                end
            end
            
            % Create weights for pareto front multiplication to turn
            % Concave fronts into convex fronts
            
            this.states.mfParetoMultipliers = this.starsAndBars(this.params.nScalarizationResolution, nObj-1);
            this.states.nScalarizations = size(this.states.mfParetoMultipliers, 1);
            
            % Initialize step sizes to half the search space
            this.states.viStep = floor((this.states.viResolution+1) / 2);
            
            this.states.viHallOfFame = [];
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % Create samples by variation along the problem axes
            miSamples = []; % Stores sampling positions
            
            % Use sample indices from hall of fame
            for sampleIndex = this.states.viHallOfFame'
                viBase = this.states.miSampleCache(sampleIndex,:);
                
                % Move base back to even grid
                viBase = round(viBase ./ this.states.viStep) .* this.states.viStep;
                
                % Confine bases to be inside the search space
                viBase = max(min(viBase, this.states.viResolution ), zeros(size(this.states.viResolution)));
                
                % Take care of initialization by adding the sample itself
                % This will be filtered out when already in cache
                miSamples = [miSamples ; viBase];
                
                % Wiggle each dimension up and down
                for iDim = [1 : numel(viBase)]
                    viSample = viBase;
                    iBase = viSample(iDim);
                    % Wiggle up and sample
                    viSample(iDim) = min(iBase + this.states.viStep(iDim), this.states.viResolution(iDim));
                    % Save to sample list
                    miSamples = [miSamples ; viSample];
                    % Wiggle down and sample
                    viSample(iDim) = max(iBase - this.states.viStep(iDim), 0);
                    % Save to sample list
                    miSamples = [miSamples ; viSample];
                end
            end
            
            % Starting condition: center of search space
            if numel(this.states.viHallOfFame) == 0
                miSamples = this.states.viStep;
            end
            
            
            % Deduplicate samples using cache
            if numel(this.states.miSampleCache)
                miDedupSamples = setdiff(miSamples, this.states.miSampleCache, 'rows');
            else
                miDedupSamples = miSamples;
            end
            
            % Filter samples and output
            samples = [];
            
            for iSample = 1:size(miDedupSamples,1)
                viSample = miDedupSamples(iSample,:);
                
                % We did not sample here before, put this sample to the cache
                this.states.miSampleCache = [this.states.miSampleCache; viSample];
                % Transform samples to floating point
                samples = [samples; viSample ./ (this.states.viResolution) .* (x_ub-x_lb) + x_lb];
            end
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            nObj = size(objectiveValues, 2);
            
            if size(objectiveValues, 1) > 0
                % Append results to objective value cache
                this.states.vfObjectiveValues = [this.states.vfObjectiveValues; objectiveValues];
                
                % Get extents of pareto front
                if numel(this.states.viHallOfFame) > 1
                    for iDim = 1:nObj
                        minObj(iDim) = min(this.states.vfObjectiveValues(:, iDim) );
                        maxObj(iDim) = max(this.states.vfObjectiveValues(this.states.viHallOfFame, iDim) );
                    end
                else
                    % Check if the bounds make sense
                    minObj = zeros(1, nObj);
                    maxObj = ones(1, nObj);
                end
                
                % Use the old fronts objective values and the currently
                % Calculated ones as the basis for decisions
                baseObjectives = [this.states.vfObjectiveValues(this.states.viHallOfFame,:); objectiveValues];
                
                % Calculate indices for base objectives
                nTotalEvals = size(this.states.vfObjectiveValues, 1);
                nLastEvals = nTotalEvals - size(objectiveValues, 1) + 1;
                baseIndices = [this.states.viHallOfFame; (nLastEvals:nTotalEvals)'];
                
                % Apply normalization to objective values
                normalizedObjectives = zeros(size(baseObjectives));
                
                for iDim = 1:size(this.states.vfObjectiveValues, 2)
                    normalizedObjectives(:, iDim) = (baseObjectives(:, iDim) - minObj(iDim))/(maxObj(iDim)-minObj(iDim));
                end
                
                viUpdateHallOfFame = [];
                this.states.viConvexFront = [];
                
                if size(normalizedObjectives, 1) <= this.params.nTrack
                    viUpdateHallOfFame = (1:size(normalizedObjectives, 1))';
                    this.states.viConvexFront = viUpdateHallOfFame;
                else
                    % Iterate over all objective dimensions
                    for iDim = 1:nObj
                        % Apply nonlinear transformations for all objective dimensions
                        for iNonlinearExponent = 1:this.states.nNonlinearizations
                            fNonlinearExponent = this.params.vfNonlinearizations(iNonlinearExponent);
                            
                            % Transform objectives by exponent
                            nonLinObjectives = normalizedObjectives;
                            nonLinObjectives(:, iDim) = sign(normalizedObjectives(:, iDim)).*(abs(normalizedObjectives(:, iDim)).^fNonlinearExponent);
                            
                            % Apply all scalarizations and apply factors to
                            % Nonlinear transformed front
                            weightedObjectives = nonLinObjectives * this.states.mfParetoMultipliers';
                            [~, sortedIndices] = sort(weightedObjectives, 1);
                            
                            % Get best indices inside the range of tracked samples
                            bestIndices = unique(sortedIndices(1:this.params.nTrack, :));
                            
                            % Append indices to hall of fame
                            viUpdateHallOfFame = union(viUpdateHallOfFame, baseIndices(bestIndices));
                            
                            % Append samples on convex front
                            this.states.viConvexFront = unique([this.states.viConvexFront; baseIndices(unique(sortedIndices(1, :)))]);
                        end
                    end
                end
                
                viUpdateHallOfFame = union(this.paretoSort(this.states.vfObjectiveValues)', viUpdateHallOfFame);
                
                % Check if we have improved the solution by comparing solutions
                bImproved = numel(setdiff(viUpdateHallOfFame, this.states.viHallOfFame)) > 0;
                this.states.viHallOfFame = viUpdateHallOfFame;
                
                if bImproved
                    terminate = false;
                    return
                end
            end
            
            % Find out which dimension has the largest step
            [iMaxStep, iMaxDim] = max(this.states.viStep);
            
            terminate = true;
            
            % If not lowest resolution reached
            if iMaxStep > 1
                terminate = false;
                
                % Reduce largest step by factor 2
                this.states.viStep(iMaxDim) = floor(this.states.viStep(iMaxDim)/ 2);
            end
            
        end
        
    end
    
    % Helper functions
    methods(Access = private)
        
        function mCombinations = starsAndBars(this, nStars, nBars)
            
            if nBars == 0
                mCombinations = nStars;
                return
            end
            
            mCombinations = [];
            
            for iStars = 0:nStars
                mSubCombination = this.starsAndBars(nStars - iStars, nBars - 1);
                mCombinations = [mCombinations; repelem(iStars, size(mSubCombination, 1), 1), mSubCombination];
            end
        end
        
    end
end