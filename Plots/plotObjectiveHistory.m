% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function fig=plotObjectiveHistory(mSamples, mResults, prob)
% PLOTOBJECTIVEHISTORY creates a plot showing objective values over
% objective function evaluations.
%
% INPUT:
% - mSamples (m by n matrix of doubles)
%   Sampled design varibles, where m is number of samples
%   and n number of design variables per sample
% - mResults (m by n matrix of doubles)
%   Sampled objective values, where m is number of samples
%   and n number of objectives per sample
% - prob (function handle)
%   Objective function
%
% OUPUT:
% - fig (figure object)
%   figure object with plot properties

res=prob([]);
if iscell(res)
    % Info about constraints
    nObj=res{1};
    nIneq=res{2};
    nEq=res{3};
else
    % No constraints
    nObj=res;
    nIneq=0;
    nEq=0;
end

fg = [];
idx_valid=ones(size(mResults, 1),1);
if ((nIneq+nEq)==0)
    % Unconstrained problem
    % Pure objective function
    fg=mResults;
else
    % Constrained problem
    for iSample = 1:size(mResults, 1)
        res = prob(mSamples(iSample, :) );
        fg(iSample, :) = res{1};
        if nIneq~=0
            if ~isempty(find(res{2}(:,:)<0))
                idx_valid(iSample)=0;
            end
        end
        if nEq~=0
            if ~isempty(find(res{3}(:,:)~=0))
                idx_valid(iSample)=0;
            end
        end
    end
end

f_valid = fg;
f_valid(idx_valid==0,1:nObj) = Inf;

% Plot convergency for each objective
fig=figure('name','Objective History');
for i =1:nObj
    subplot(1, nObj, i);
    hold on
    % Pure objective function values
    plot([1:size(mResults, 1)],fg(:,i), '.k')
    % Current minimal valid objective funcion value 
    plot([1:size(mResults, 1)],cummin(f_valid(:,i)), '-r','LineWidth',2)
%     % Current minimal extended objective funcion value 
%     plot([1:size(mResults, 1)],cummin(mResults(:,i)), '-r','LineWidth',2)
    xlabel('Objective Function Evaluations')
    ylabel(['Objective Value ', num2str(i)]);
    
    legend('Samples', 'Best Objective Value', 'Location', 'northeast');
end

end