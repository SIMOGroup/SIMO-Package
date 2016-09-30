function PlotKnts(NURBS)
% function PlotKnts(NURBS)

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

if NURBS.Dim == 1
    EvalKnts = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, NURBS.uqKntVect);
    Weights = reshape(EvalKnts(4, :), 1, []);
    Pts = bsxfun(@rdivide, EvalKnts, Weights);
    plot3(Pts(1, :), Pts(2, :),  Pts(3, :),'s','MarkerEdgeColor', 'g',...
        'MarkerFaceColor', 'g', 'MarkerSize', 5);
elseif NURBS.Dim == 2
    ParaPts{1} = linspace(0, 1, 101); % parametric points
    ParaPts{2} = linspace(0, 1, 101); % parametric points
    % Plot the knots
    Curvs{1} = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {NURBS.uqKntVect{1}, ParaPts{2}});
    Curvs{2} = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {ParaPts{1}, NURBS.uqKntVect{2}});
    
    for i = 1 : numel(NURBS.uqKntVect{1})
        EtaCurv = reshape(Curvs{1}(1 : 3, i, :), 3, []);
        Weights = reshape(Curvs{1}(4, i, :), 1, []);
        EtaCurv = bsxfun(@rdivide, EtaCurv, Weights);
        plot3(EtaCurv(1, :), EtaCurv(2, :), EtaCurv(3, :),...
            'linewidth', 1.5, 'color', [95 137 60] ./ 250);
        %         plot3(EtaCurv(1, :), EtaCurv(2, :), EtaCurv(3, :),...
        %             'linewidth', 1.5, 'color', 'k');
    end
    for j = 1 : numel(NURBS.uqKntVect{2})
        XiCurv = reshape(Curvs{2}(1 : 3, :, j), 3, []);
        Weights = reshape(Curvs{2}(4, :, j), 1, []);
        XiCurv = bsxfun(@rdivide, XiCurv, Weights);
        plot3(XiCurv(1, :), XiCurv(2, :), XiCurv(3, :),...
            'linewidth', 1.5, 'color', [95 137 60] ./ 250);
        %         plot3(XiCurv(1, :), XiCurv(2, :), XiCurv(3, :),...
        %             'linewidth', 1.5, 'color', 'k');
    end
elseif NURBS.Dim == 3
    for dir = 1 : 3
        for side = [0 1]
            Surf = NURBSExtract(NURBS, dir, side);
            PlotKnts(Surf);
        end
    end
end
end