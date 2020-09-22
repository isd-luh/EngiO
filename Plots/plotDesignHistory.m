% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function fig=plotDesignHistory(mSamples)
% PLOTDESIGNHISTORY creates a plot showing design variable values over
% objective function evaluations.
%
% INPUT:
% - mSamples (m by n matrix of doubles)
%   Sampled design varibles, where m is number of samples
%   and n number of design variables per sample
%
% OUPUT:
% - fig (figure object)
%   figure object with plot properties

% Plot convergency of each design variable
fig=figure('name','Design History');
hold on
for i =1:size(mSamples, 2)
    subplot(1, size(mSamples, 2), i)
    plot(mSamples(:, i), '.k')
    xlabel('Objective Function Evaluations')
    ylabel(['Design Variable ', num2str(i)])
end

end