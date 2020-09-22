% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function fig=plotConstraintHistory(mSamples, mResults, prob)
% PLOTCONSTRAINTHISTORY creates a plot showing constraint equation values over
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
    error('No constraints');
end

g = [];
h = [];
for iSample = 1:size(mSamples, 1)
    res = prob(mSamples(iSample, :) );
    if nIneq~=0
        g(iSample, :) = res{2};
    end
    if nEq~=0
        h(iSample, :) = res{3};
    end
end

fig=figure('name','Constraint History');
for i =1:nIneq
    subplot(1, nIneq+nEq, i);
    plot(1:size(mResults,1), g(:,i),'.k')
    xlabel('Objective Function Evaluations')
    ylabel(['Constraint Value ',num2str(i)]);
end
for i =1:nEq
    subplot(1, nIneq+nEq, nIneq+i);
    plot(1:size(mResults,1), h(:,i),'.k')
    xlabel('Objective Function Evaluations')
    ylabel(['Constraint Value ',num2str(nIneq+i)]);
end

end