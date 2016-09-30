function D = getElastMat(E, nu, lab)

% Evaluate elasticity matrix

%{
Copyright (C) <2014-2016>  <Khanh Chau-Nguyen, Hung Nguyen-Xuan>

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

if strcmp(lab, 'PlaneStress')
    % Plane stress state
    D = E / (1 - nu ^ 2) * [1 nu 0;
        nu 1 0;
        0   0  (1 - nu) / 2];
elseif strcmp(lab, 'PlaneStrain')
    D = (E / ((1 + nu) * (1 - 2 * nu))) * [1 - nu, nu, 0;
        nu, 1 - nu, 0;
        0, 0, (1 - 2 * nu) / 2];
elseif strcmp(lab, '3D') % 3D solid
    D = zeros(6, 6);
    D(1 : 3, 1 : 3) = E / (1 + nu)/(1 - 2 * nu) * [1 - nu     nu     nu;
        nu 1 - nu     nu;
        nu     nu 1 - nu];
    D(4 : 6, 4 : 6) = E / (2 * (1 + nu)) * eye(3);
else
    error('"lab" must be "PlaneStress", "PlaneStrain", or "3D"');
end
end