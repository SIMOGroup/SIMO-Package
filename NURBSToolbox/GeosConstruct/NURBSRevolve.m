function surf = NURBSRevolve(curve,pnt,vec,theta)

%
% NRBREVOLVE: Construct a NURBS surface by revolving a NURBS curve, or
%  construct a NURBS volume by revolving a NURBS surface.
%
% Calling Sequence:
%
%   srf = nrbrevolve(crv,pnt,vec[,ang])
%
% INPUT:
%
%   crv		: NURBS curve or surface to revolve, see nrbmak.
%
%   pnt		: Coordinates of the point used to define the axis
%               of rotation.
%
%   vec		: Vector defining the direction of the rotation axis.
%
%   ang		: Angle to revolve the curve, default 2*pi
%
% OUTPUT:
%
%   srf		: constructed surface or volume
%
% Description:
%
%   Construct a NURBS surface by revolving the profile NURBS curve around
%   an axis defined by a point and vector.
%
% Examples:
%
%   Construct a sphere by rotating a semicircle around a x-axis.
%
%   crv = nrbcirc(1.0,[0 0 0],0,pi);
%   srf = nrbrevolve(crv,[0 0 0],[1 0 0]);
%   nrbplot(srf,[20 20]);
%
% NOTE:
%
%   The algorithm:
%
%     1) vectrans the point to the origin (0,0,0)
%     2) rotate the vector into alignment with the z-axis
%
%     for each control point along the curve
%
%     3) determine the radius and angle of control
%        point to the z-axis
%     4) construct a circular arc in the x-y plane with
%        this radius and start angle and sweep angle theta
%     5) combine the arc and profile, coefs and weights.
%
%     next control point
%
%     6) rotate and vectrans the surface back into position
%        by reversing 1 and 2.
%
%
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2010 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if (nargin < 3)
    error('Not enough arguments to construct revolved surface');
end

if (nargin < 4)
    theta = 2.0*pi;
end

if (iscell (curve.KntVect) && numel(curve.KntVect) == 3)
    error('The function nrbrevolve is not yet ready to create volumes')
end

% Translate curve the center point to the origin
if isempty(pnt)
    pnt = zeros(3,1);
end

if length(pnt) ~= 3
    error('All point and vector coordinates must be 3D');
end

% Translate and rotate the original curve or surface into alignment with the z-axis
T  = vectrans(-pnt);
angx = vecangle(vec(1),vec(3));
RY = vecroty(-angx);
vectmp = RY*[vecnorm(vec(:));1.0];
angy = vecangle(vectmp(2),vectmp(3));
RX = vecrotx(angy);
curve = NURBSTForm(curve,RX*RY*T);

% Construct an arc
arc = NURBSCirc(1.0,[],0.0,theta);

if (numel(curve.KntVect) > 1)
    % Construct the revolved volume
    CtrlPts = zeros([4 arc.NCtrlPts curve.NCtrlPts]);
    angle = squeeze (vecangle(curve.CtrlPts4D(2,:,:),curve.CtrlPts4D(1,:,:)));
    radius = squeeze (vecmag(curve.CtrlPts4D(1:2,:,:)));
    for i = 1:curve.NCtrlPts(1)
        for j = 1:curve.NCtrlPts(2)
            CtrlPts(:,:,i,j) = vecrotz(angle(i,j))*vectrans([0.0 0.0 curve.CtrlPts4D(3,i,j)])*...
                vecscale([radius(i,j) radius(i,j)])*arc.CtrlPts4D;
            CtrlPts(4,:,i,j) = CtrlPts(4,:,i,j)*curve.CtrlPts4D(4,i,j);
        end
    end
    surf = CreateNURBS([arc.KntVect, curve.KntVect], CtrlPts);
else
    % Construct the revolved surface
    CtrlPts = zeros(4, arc.NCtrlPts, curve.NCtrlPts);
    angle = vecangle(curve.CtrlPts4D(2,:),curve.CtrlPts4D(1,:));
    radius = vecmag(curve.CtrlPts4D(1:2,:));
    for i = 1 : curve.NCtrlPts
        CtrlPts(:,:,i) = vecrotz(angle(i))*vectrans([0.0 0.0 curve.CtrlPts4D(3,i)])*...
            vecscale([radius(i) radius(i)])*arc.CtrlPts4D;
        CtrlPts(4,:,i) = CtrlPts(4,:,i)*curve.CtrlPts4D(4,i);
    end
    surf = CreateNURBS([arc.KntVect, curve.KntVect], CtrlPts);
end

% Rotate and vectrans the surface back into position
T = vectrans(pnt);
RX = vecrotx(-angy);
RY = vecroty(angx);
surf = NURBSTForm(surf,T*RY*RX);

end

