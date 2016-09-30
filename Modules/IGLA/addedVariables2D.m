function [A2] = addedVariables2D(sdof, k1, k2, Bstrain, totalGP, stressState)       % Plane Strain case)

% function [A2] = addedVariables2D(sdof, k1, k2, Bstrain, totalGP, stressState)
%{
Copyright (C) <2014-2016>  <Hung Nguyen-Xuan, Khanh Chau-Nguyen>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

A2 = sparse(k2, sdof);
count = 1;
if ( strcmp(stressState,'PlaneStress') )       % Plane Strain case
    for i = 1 : totalGP
        Bx = Bstrain{i}(1, :); % see how to store this value in evaluated shape function
        By = Bstrain{i}(2, :);
        Bxy = Bstrain{i}(3, :);
        A2(count,:) = 2/sqrt(3)*Bx + 1/sqrt(3)*By;
        A2(count + 1, :) = By;
        A2(count + 2, :) = 1 / sqrt(3)*Bxy;
        count = count + 3;
        clear Bx By Bxy;
    end
else
    for i = 1 : totalGP
        Bx = Bstrain{i}(1, :); % see how to store this value in evaluated shape function
        By = Bstrain{i}(2, :);
        Bxy = Bstrain{i}(3, :);
        A2(count, :) = Bx - By;
        A2(count + 1, :) = Bxy;
        count = count + 2;
        clear Bx By Bxy;
    end
end
A2 = [A2, sparse(k2, k1), -sparse(1 : k2, 1 : k2, ones(k2, 1))];
% A2 = sparse([[A2], [zeros(k2,k1)], [-eye(k2)]]);