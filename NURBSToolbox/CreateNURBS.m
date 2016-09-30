function NURBS = CreateNURBS(KntVect, CtrlPts)
% NURBS = CreateNURBS(KntVects, CtrlPts)
% -------------------------------------------------------------------------
% Gathering and constructing data for a NURBS patch (1D, 2D, 3D)
%--------------------------------------------------------------------------
% Input:
%       KntVect: knot vectors stored in cell format
%       CtrlPts: control points stored in homogenous co-ordinates, i.e. 
% location of an arbitrary control point is identified by [x1*w1; y1*w1; z1*w1; w1]
%--------------------------------------------------------------------------
% Output:
%       NURBS: a NURBS structure, including
%           NURBS.NDS       : number of dimensional space
%           NURBS.KntVect   : knot vector(s) stored in cell format
%           NURBS.uqKntVect : unique knot values of knot vector(s) in cell format
%           NURBS.KntMult   : multiplicities of knot values
%           NURBS.CtrlPts4D : control points co-ordinates in 4D space
%           NURBS.CtrlPts3D : control points co-ordinates projected into 3D space
%           NURBS.Weights   : weights of control points
%           NURBS.Dim       : number of dimensions of the NURBS patch
%           NURBS.NCtrlPts  : number of control points in each direction
%           NURBS.Order     : degree of basis functions in each direction
%           NURBS.NNP       : total number of control points, "NP" is an abbreviation for nodal points
%--------------------------------------------------------------------------

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

assert(iscell(KntVect), 'Input knots must be in cell format');
assert(numel(KntVect) >= 1);
assert(numel(KntVect) <= 3);
Dim = numel(KntVect); % The dimensionality of the model
for i = 1 : Dim
    assert(size(KntVect{i}, 1) == 1)
    assert(numel(KntVect{i}) >= 4, ...
        'Number of knot values must be equal or greater than 4')
end
% ------------------------------------
W = CtrlPts(4, :, :, :);
CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :, :, :), W);

NPts = size(CtrlPts);
NCtrlPts = zeros(1, Dim);
p = zeros(1, Dim);
uqKntVect = cell(1, Dim);
uqKntsIdcs = cell(1, Dim);
KntMult = cell(1, Dim);
for i = 1 : Dim
    NCtrlPts(i) = NPts(i + 1);
    p(i) = numel(KntVect{i}) - NCtrlPts(i) - 1;
    assert(p(i) > 0, 'Degree of the spline is negative, please check input parameters!');
    uqKntsIdcs{i} = [true, diff(KntVect{i}) > 0];
    uqKntVect{i}  = KntVect{i}(uqKntsIdcs{i});
    KntMult{i} =  diff([find(uqKntsIdcs{i} == 1), numel(KntVect{i}) + 1]);
end
% check number of dimensional space
if all(abs(CtrlPts(2 : 3, :, :, :)) <= eps)
    NURBS.NDS = 1;
elseif all(abs(CtrlPts(3, :, :, :)) <= eps)
    NURBS.NDS = 2;
else
    NURBS.NDS = 3;
end
NURBS.KntVect = KntVect;
NURBS.uqKntVect = uqKntVect; % unique knot values
NURBS.KntMult = KntMult; % multiplicities of knot values
NURBS.CtrlPts4D = CtrlPts; % control points in 4D space
NURBS.CtrlPts3D = CtrlPts3D; % control points projected into 3D space
NURBS.Weights = W;
NURBS.Dim = Dim;
% number of control points in each direction
NURBS.NCtrlPts = NCtrlPts;
NURBS.Order = p;
% number of total control points
NURBS.NNP = prod(NCtrlPts); %"NP" is an abbreviation for nodal points
end