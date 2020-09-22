% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

classdef Grid < Optimizer
    % GRID is a grid sampling implementation derived from Optimizer.
    %
    
    methods
        function this = Grid(states_init, options)
            
            % there are no parameters
            this.defaultparams = struct();
            
        end
        
        function initialize(this, nObj, vartype, x_lb, x_ub)
            
            % calculate resolution based on the maximum number of
            % evaluations allowed
            resolution = ceil(this.options.maxEvals^(1/numel(x_lb) ));
            
            % Initialize all states
            this.states = struct('nResolution', resolution);
            this.states.vBest = x_lb;
            
        end
        
        function samples = generateSamples(this, nObj, vartype, x_lb, x_ub)
            
            % get number of dimensions
            nDims = numel(x_lb);
            
            % Generate grid samples
            mGrid = 0 : double(this.states.nResolution-1);
            for iDim = [1:nDims-1]
                mGrid = [repelem(mGrid(1,:), this.states.nResolution);repmat(mGrid, 1, this.states.nResolution)];
            end
            nCurrSamples = numel(mGrid(1,:));
            
            % Squash grid to boundaries
            mNorm = diag(x_ub - x_lb) / (double(this.states.nResolution)-1);
            mOffset = repmat(x_lb, nCurrSamples, 1);
            samples = mGrid' * mNorm + mOffset;
            
            % Integer variables are rounded
            for j = 1:size(samples, 1)
                if ~isempty(vartype)
                    for m = 1:length(x_lb)
                        if strcmpi(vartype{m}, 'int')
                            samples(j, m) = round(samples(j, m));
                        end
                    end
                end
            end
            
        end
        
        function terminate = processResults(this, samples, objectiveValues)
            
            % terminate optimization run after first iteration
            terminate = true;
            
        end
        
    end
end