% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef GlobalPatternDepthFirst < Optimizer
    % GLOBALPATTERNDEPTHFIRST is a single-objective deterministic algorithm
    % with depth-first search instead of breadth first. It is
    % derived from Optimizer.
    %
    % LITERATURE:
    %   B. Hofmeister, M. Bruns, R. Rolfes (2019) Finite element model updating
    %   using deterministic optimisation: A global pattern search approach,
    %   Engineering Structures (195)
    methods
        function this = GlobalPatternDepthFirst()
            this.defaultparams = struct('nTrack', 15, ...   % local minima to track
                'nResolution', 2^20 ... % resolution of the grid
                );
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            nDims = numel(x_lb);
            
            this.states = struct('viResolution', zeros(nDims), ...  % Allocate space for resolution vector
                'iTrack', 1, ...                   % Number of currently tracked optima
                'miHallOfFame', [], ...
                'vfHallOfFame', [], ...
                'miActualSamples', [], ...
                'bCachedMin', false, ...
                'miSampleCache', [], ...            % Cache for sample coordinates
                'vfObjectiveValues', []);           % Cache for objective values
            
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
            
            % Initialize step sizes to half the search space
            this.states.viStep = floor((this.states.viResolution+1) / 2);
            
            this.states.miHallOfFame = repmat(this.states.viStep, this.states.iTrack, 1);
            this.states.vfHallOfFame = repelem(Inf, this.states.iTrack);
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % First the algorithm chooses the bases for a variation
            
            % Initialize and use all previous locations initially
            
            % Now the bases are fixed and the samples are created by
            % Variation along the problem axes
            
            miSamples = []; % stores sampling positions
            
            for viBaseT = this.states.miHallOfFame'
                viBase = viBaseT';
                
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
            
            % Filter samples and output
            samples = [];
            
            oldCache = this.states.miSampleCache;
            this.states.miActualSamples = [];
            this.states.bCachedMin = false;
            
            for viSampleT = miSamples'
                viSample = viSampleT';
                % Check if we already sampled this point
                if numel(oldCache) > 1
                    [mem, index] = ismember(viSample, oldCache, 'rows');
                    if mem
                        % Modified hall of fame update based on old cache entries for depth first
                        if ~ismember(viSample, this.states.miHallOfFame, 'rows');
                            [maxval, maxIndices] = max(this.states.vfHallOfFame);
                            
                            if this.states.vfObjectiveValues(index) < maxval
                                this.states.vfHallOfFame(maxIndices(1)) = this.states.vfObjectiveValues(index);
                                this.states.miHallOfFame(maxIndices(1),:) = viSample;
                                this.states.bCachedMin = true;
                            end
                        end
                        continue
                    end
                    
                    
                end
                if numel(this.states.miActualSamples) > 1
                    if ismember(viSample, this.states.miActualSamples, 'rows')
                        continue
                    end
                end
                % We did not sample here before, put this sample to the cache
                this.states.miSampleCache = [this.states.miSampleCache; viSample];
                this.states.miActualSamples = [this.states.miActualSamples; viSample];
                % Transform samples to floating point
                samples = [samples; viSample ./ (this.states.viResolution) .* (x_ub-x_lb) + x_lb];
            end
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % Append results to objective value cache
            this.states.vfObjectiveValues = [this.states.vfObjectiveValues; objectiveValues];
            
            for iResult = 1:numel(objectiveValues)
                
                if ismember(this.states.miActualSamples(iResult,:), this.states.miHallOfFame, 'rows') && ~any(this.states.vfHallOfFame == Inf)
                    continue
                end
                
                [maxval, maxIndices] = max(this.states.vfHallOfFame);
                
                if objectiveValues(iResult) < maxval
                    this.states.vfHallOfFame(maxIndices(1)) = objectiveValues(iResult);
                    this.states.miHallOfFame(maxIndices(1),:) = this.states.miActualSamples(iResult,:);
                    terminate = false;
                end
            end
            
            if this.states.bCachedMin
                terminate = false;
            end
            
            % Find out which dimension has the largest step
            [iMaxStep, iMaxDim] = max(this.states.viStep);
            
            terminate = true;
            
            % If not lowest resolution reached
            if iMaxStep > 1
                terminate = false;
                
                % Reduce largest step by factor 2
                this.states.viStep(iMaxDim) = floor(this.states.viStep(iMaxDim)/ 2);
            elseif this.states.iTrack <= this.params.nTrack
                % Increase number of tracked optima
                this.states.iTrack = this.states.iTrack + 1;
                
                % Initialize step sizes to half the search space
                this.states.viStep = floor((this.states.viResolution+1) / 2);
                
                this.states.miHallOfFame = repmat(this.states.viStep, this.states.iTrack, 1);
                this.states.vfHallOfFame = repelem(Inf, this.states.iTrack);
                
                terminate = false;
            end
            
        end
        
    end
end