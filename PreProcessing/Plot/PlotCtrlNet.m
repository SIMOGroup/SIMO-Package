function PlotCtrlNet(NURBS)
% function PlotCtrlNet(NURBS)

if NURBS.Dim == 1
    CtrlPts3D = NURBS.CtrlPts3D(1 : 3, :);
    plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
elseif NURBS.Dim == 2
    % Plot control points along y-direction
    for i = 1 : NURBS.NCtrlPts(1)
        CtrlPts3D = reshape(NURBS.CtrlPts3D(:, i, :), 3, []);
        plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--');
    end
    
    % Plot control points along x-direction
    for j = 1 : NURBS.NCtrlPts(2)
        CtrlPts3D = reshape(NURBS.CtrlPts3D(:, :, j), 3, []);
        plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
    end
elseif NURBS.Dim == 3
    for i = 1 : NURBS.NCtrlPts(1)
        for j = 1 : NURBS.NCtrlPts(2)
            CtrlPts3D = reshape(NURBS.CtrlPts3D(:, i, j, :), 3, []);
            plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
        end
        for k = 1 : NURBS.NCtrlPts(3)
            CtrlPts3D = reshape(NURBS.CtrlPts3D(:, i, :, k), 3, []);
            plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
        end
    end
    for j = 1 : NURBS.NCtrlPts(2)
        for k = 1 : NURBS.NCtrlPts(3)
            CtrlPts3D = reshape(NURBS.CtrlPts3D(:, :, j, k), 3, []);
            plot3(CtrlPts3D(1, :), CtrlPts3D(2, :), CtrlPts3D(3, :),'k--')
        end
    end    
end
end