function [L2Norm,H1Norn,ENorm] = calErrorNorm2D(Mesh, NURBS, E, U, ExactSolu)
% function [DNorm, ENorm] = calErrorNorm(Mesh, NURBS, E, U, ExactSolu)
% ----------------------------------------------------------------------
% Calculate displacement norm and energy norm for bar problem
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure of the computational domain
%       NURBS: NURBS structure
%       E: Young's modulus
%       A: area of cross section
%       U: nodal displacement vector
% ---------------------------------------------------------------------
% Output:
%       DNorm: displacement error norm
%       ENorm: energy error norm
%
%       L2 norm: 
%       H1 Norm;
%       
% --------------------------------------------------------------------

NGPs = NURBS.Order + 1; % number of gauss points

[Jx, Wx, ~, Nx] = calcDersBasisFunsAtGPs(NURBS.Order(1), NURBS.NCtrlPts(1), NURBS.KntVect{1}, 2, NGPs(1), Mesh.NElDir(1));
[Jy, Wy, ~, Ny] = calcDersBasisFunsAtGPs(NURBS.Order(2), NURBS.NCtrlPts(2), NURBS.KntVect{2}, 2, NGPs(2), Mesh.NElDir(2));
N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);
N2 = zeros(3, Mesh.NEN);

Weights = NURBS.Weights;
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

Weights = reshape(NURBS.Weights, 1, []);
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

L2Norm_l = 0;
H1Norn_l = 0;
ENorm_l = 0;

dispNorm = 0;
energyNorm = 0;

% Loop over elements (knot spans)
for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        for qy = 1 : NGPs(2)
            for qx = 1 : NGPs(1)
                k = 1;
                for j = 1 : NURBS.Order(2) + 1
                    for i = 1 : NURBS.Order(1) + 1
                        N0(k)    = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1); %N.N
                        N1(1, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1); %dN/dxi.N 
                        N1(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2); %N.dN/deta
                        N2(1, k) = Nx(ex, qx, i, 3) * Ny(ey, qy, j, 1); % d2N/dxi2
                        N2(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 3); % d2N/deta2
                        N2(3, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 2); % d2N/(dxi.deta)
                        k = k + 1;
                    end
                end
                [R0, R1, R2] = Rationalize(Weights(Mesh.El(e, :)), N0, N1, N2);
                J2 = Jx(ex) * Jy(ey);
                W = Wx(qx) * Wy(qy);
                % Gradient of mapping from parameter space to physical space
                jacob  = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                jacob2 = R2 * CtrlPts(Mesh.El(e, :), 1 : 2);
                
                J1 = abs(det(jacob));
                
                dxdxi = jacob(1, 1); 
                dydxi = jacob(1, 2);
                dxdeta = jacob(2, 1); 
                dydeta = jacob(2, 2);
                
                j33 = [dxdxi^2, dydxi^2, 2*dxdxi*dydxi;
                    dxdeta^2, dydeta^2, 2*dxdeta*dydeta;
                    dxdxi*dxdeta, dydxi*dydeta, dxdxi*dydeta + dxdeta*dydxi];
            
                % Jacobian inverse and spatial 1st and 2nd derivatives
                % Compute derivatives of basis functions w.r.t
                % physical coordinates
                
                dRdx = jacob \ R1;
                dR2dx = j33 \ (R2 - jacob2*dRdx);
                
                % physical coords of gauss points
                Pts = R0 * CtrlPts(Mesh.El(e, :), :);         
                x = Pts(1);   
                y = Pts(2);   
                %z = Pts(1);  
                % Numerical solution
                strain = dRdx * U(Mesh.El(e, :));                        
                duh = E(x) * strain;
                uh = R0 * U(Mesh.El(e, :));
                % Exact solution
                r = sqrt(x^2+y^2);
                ue = ur.disp(r);
                due = ur.strain(r);
                %Erorr calculation
                L2Norm = L2Norm + (ue - uh)'*(ue - uh)*J1*J2*W;     
                %H1Norn = ````;              %
                ENorm = ENorm + (due - duh)'*(due - duh)*J1*J2*W;
            end
        end   
    end
end

L2Norm = sqrt(L2Norm);
H1Norm = sqrt(H1Norm);
ENorm = sqrt(energyNorm);
end
