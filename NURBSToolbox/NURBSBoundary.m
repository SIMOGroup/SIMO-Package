function Bdry = NURBSBoundary(NURBS, iSide)
% function Bdry = NURBSBoundary(NURBS, iSide)
% --------------------------------------------------------------
% Extract the boundary of a NURBS structure
% --------------------------------------------------------------
% For 1D case:
% 
%         iSide = 1 o------------o iSide = 2     +----> u
%          (u = 0)                  (u = 1)
% --------------------------------------------------------------
% For 2D case:
%                 iSide = 4
%                  (v = 1)
%              o------------o                   v
%              |            |                   ^
%              |            |                   |
%   iSide = 1  |            | iSide = 2         |
%    (u = 0)   |            |  (u = 1)          +----> u
%              |            |
%              o------------o
%                 iSide = 3
%                  (v = 0)
% --------------------------------------------------------------
% For 3D case:
% 
%                    o------------o
%                   /| iSide = 6 /|
%                  / | (w = 1)  / |          w
%                 o------------o  iSide = 2  ^  v
%                 |  iSide = 4 |  |(u = 1)   | /
%       iSide = 1 |  | (v = 1) |  |          |/
%        (u = 0)  |  o-------- | -o          +----> u
%                 | /  iSide = 3 / 
%                 |/    (v = 0)|/
%                 o------------o
%                    iSide = 5   
%                     (w = 0)
% --------------------------------------------------------------
    
Idx = floor((iSide + 1) / 2);
if mod(iSide, 2) == 0
    Bdry = NURBSExtract(NURBS, Idx, 1);
else
    Bdry = NURBSExtract(NURBS, Idx, 0);
end
end