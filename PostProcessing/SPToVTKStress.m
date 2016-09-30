function SPToVTKStress(NURBS, d, ParaPts, E, nu, lab, filename)
% SPToVTKStress(NURBS, d, ParaPts, E, nu, lab, filename)
% -------------------------------------------------------------------------
% Export stresses to *.vtk
% -------------------------------------------------------------------------
% Input:
%       NURBS: NURBS structure
%       d: displacement vector
%       ParaPts: parameter points
%       E: Young's modulus
%       nu: Poisson's ratio
%       lab: label to identify stress state 
%           (PlaneStress, PlaneStrain, or 3D)
%       filename: name of vtk file
% -------------------------------------------------------------------------
% Output:
%       filename.vtk
% -------------------------------------------------------------------------

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

D = getElastMat(E, nu, lab);
F1 = reshape(d, NURBS.NNP, [])';
F2 = reshape(reshape(d, NURBS.NNP, [])', [size(F1, 1), NURBS.NCtrlPts]);
g = GradEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts, F2);
NPts = cellfun(@numel, ParaPts);
Cw = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, ParaPts);
if NURBS.Dim == 2
    stress = zeros([3, NPts]);
    for i = 1 : NPts(1)
        for j = 1 : NPts(2)
            strain(1) = g(1, 1, i, j); % \epsilon_{xx}
            strain(2) = g(2, 2, i, j); % \epsilon_{yy}
            strain(3) = g(2, 1, i, j) + g(1, 2, i, j); % 2\epsilon_{xy}
            stress(:, i, j) = D * strain';
        end
    end
    S.xx = stress(1, :, :);
    S.yy = stress(2, :, :);
    if strcmp(lab, 'PlaneStrain')
        S.zz = nu * (S.xx + S.yy);
    elseif strcmp(lab, 'PlaneStress')
        S.zz = zeros(size(S.xx));
    else
        error('"lab" must be "PlaneStrain" or "PlaneStress"');
    end
    S.xy = stress(3, :, :);
    S.yz = zeros(size(S.xx));
    S.xz = zeros(size(S.xx));
    
    Weights = Cw(4, :, :);
    C = bsxfun(@rdivide, Cw(1 : 3, :, :), Weights);
else
    stress = zeros([6, NPts]);
    for i = 1 : NPts(1)
        for j = 1 : NPts(2)
            for k = 1 : NPts(3)
                strain(1) = g(1, 1, i, j, k); % \epsilon_{xx}
                strain(2) = g(2, 2, i, j, k); % \epsilon_{yy}
                strain(3) = g(3, 3, i, j, k); % \epsilon_{zz}
                strain(4) = g(2, 1, i, j, k) + g(1, 2, i, j, k); % 2\epsilon_{xy}
                strain(5) = g(3, 2, i, j, k) + g(2, 3, i, j, k); % 2\epsilon_{yz}
                strain(6) = g(3, 1, i, j, k) + g(1, 3, i, j, k); % 2\epsilon_{xz}
                stress(:, i, j, k) = D * strain';
            end
        end
    end
    S.xx = stress(1, :, :, :);
    S.yy = stress(2, :, :, :);
    S.zz = stress(3, :, :, :);
    S.xy = stress(4, :, :, :);
    S.yz = stress(5, :, :, :);
    S.xz = stress(6, :, :, :);
    
    Weights = Cw(4, :, :, :);
    C = bsxfun(@rdivide, Cw(1 : 3, :, :, :), Weights);
end
% Calculate vonMises stress
firstTerm = 1 / 2 * ((S.xx - S.yy) .^ 2 + ...
    (S.yy - S.zz) .^ 2 + ...
    (S.zz - S.xx) .^ 2);
secondTerm = 3 * (S.xy .^ 2 + S.yz .^ 2 + S.xz .^ 2);
vonMises = sqrt(firstTerm + secondTerm);

str1 = cat (2,'<?xml version="1.0"?> \n', ...
    '<VTKFile type="StructuredGrid" version="0.1"> \n', ...
    '<StructuredGrid WholeExtent="0 %d 0 %d 0 %d"> \n', ...
    '<Piece Extent="0 %d 0 %d 0 %d"> \n', ...
    '<PointData Scalars="Stress components">\n');

str2 = ...
    '<DataArray type="Float32" Name="%s" format="ascii" NumberOfComponents="1"> \n';

str3 = '</DataArray> \n';

str4 = cat (2,'</PointData> \n', ...
    '<Points> \n', ...
    '<DataArray type="Float32" NumberOfComponents="3"> \n');

str5 = cat (2, '\n', ...
    '</DataArray>\n', ...
    '</Points> \n', ...
    '</Piece> \n', ...
    '</StructuredGrid> \n', ...
    '</VTKFile> \n');

if (NURBS.Dim == 2)
    NPts(3) = 1;
    C(3, :, :) = 0;
end

if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.vts'))
    filename = cat (2, filename, '.vts');
end

fid = fopen (filename, 'w');
if (fid < 0)
    error ('SPToVTKStress: could not open file %s', filename);
end

fprintf (fid, str1, ...
    NPts(1) - 1, NPts(2) - 1, NPts(3) - 1, ...
    NPts(1) - 1, NPts(2) - 1, NPts(3) - 1);

%--------------------------------------------------------------------------
fprintf (fid, str2, 'sigma_xx');
fprintf (fid, '%g ', reshape(squeeze(S.xx), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str2, 'sigma_yy');
fprintf (fid, '%g ', reshape(squeeze(S.yy), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str2, 'sigma_zz');
fprintf (fid, '%g ', reshape(squeeze(S.zz), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str2, 'sigma_xy');
fprintf (fid, '%g ', reshape(squeeze(S.xy), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str2, 'sigma_yz');
fprintf (fid, '%g ', reshape(squeeze(S.yz), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str2, 'sigma_xz');
fprintf (fid, '%g ', reshape(squeeze(S.xz), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str2, 'vonMises');
fprintf (fid, '%g ', reshape(squeeze(vonMises), [], 1));
fprintf (fid, str3);
%
fprintf (fid, str4);
fprintf (fid, '%g ', C(:));
fprintf (fid, str5);
fclose (fid);
end
