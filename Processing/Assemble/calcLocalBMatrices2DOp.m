function [Bstrain, WeJ, countGP]= calcLocalBMatrices2DOp(Mesh, NURBS, stressState)
% function [Bstrain, WeJ, countGP] = calcLocalBMatrices2DOp(Mesh, NURBS)
% -------------------------------------------------------------------
% Calculate two dimensional local B matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
% -------------------------------------------------------------------
% Output:
%       Bstrain: matrix stores local B matrices at all Gauss points
%       (size(Bstrain) = NEl * NGPs
%       WeJ (Nel*NGPs): vector stores J1*J2*W through all Gauss points
%       countGP: total numeber of Gauss points over entire problem
% -------------------------------------------------------------------

% % config for p = q = 2, regularity = 1
% [upara, Wx] = quadrule(2*NURBS.Order(1), 0, Mesh.NElDir(1), 0, 1);
% [vpara, Wy] = quadrule(2*NURBS.Order(2), 0, Mesh.NElDir(2), 0, 1);

% config for p = q = 3, regularity = 2
[upara, Wx] = quadrule(2*NURBS.Order(1), 1, Mesh.NElDir(1), 0, 1);
[vpara, Wy] = quadrule(2*NURBS.Order(2), 1, Mesh.NElDir(2), 0, 1);

% corners = NURBSEval(NURBS, {upara, vpara});
% x = reshape(corners(1, :, :), 1, [])';
% y = reshape(corners(2, :, :), 1, [])';
% plot(x, y, 'kx'); hold on;

spanIdxX = FindSpan(NURBS.NCtrlPts(1), NURBS.Order(1), upara, NURBS.KntVect{1});
spanIdxY = FindSpan(NURBS.NCtrlPts(2), NURBS.Order(2), vpara, NURBS.KntVect{2});

Nx = DersBasisFuns(spanIdxX, upara, NURBS.Order(1), 1, NURBS.KntVect{1});
Ny = DersBasisFuns(spanIdxY, vpara, NURBS.Order(2), 1, NURBS.KntVect{2});

% call coordinates of isoparametric 4-node element (parent Q4 element)
Q4=[-1 -1;1 -1;1 1;-1 1];

N0 = zeros(1, Mesh.NEN);
N1 = zeros(NURBS.Dim, Mesh.NEN);

Weights = reshape(NURBS.Weights, 1, []);
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

countGP = 0;   % count the number of Gauss points
Bstrain = cell(numel(upara) * numel(vpara), 1);
WeJ = zeros(numel(upara) * numel(vpara), 1);

tmp1 = [true; diff(spanIdxX) > 0];
tmp2 = [true; diff(spanIdxY) > 0];
spanMultU = diff([find(tmp1 == 1); numel(spanIdxX) + 1]);
spanMultV = diff([find(tmp2 == 1); numel(spanIdxY) + 1]);
idcsQPU = cumsum([1; spanMultU]);
idcsQPV = cumsum([1; spanMultV]);

%if ( strcmp(stressState,'PlaneStress') )       % Plane Stress case
    for e = 1 : Mesh.NEl
        [ex, ey] = ind2sub(Mesh.NElDir, e);
        
        getDofNEl(1 : Mesh.NEN) = Mesh.El(e, :);
        getDofNEl(Mesh.NEN + 1 : 2*Mesh.NEN) = Mesh.El(e, :) + Mesh.NDof / Mesh.Dof;
        
        for qy = idcsQPV(ey) : idcsQPV(ey + 1) - 1
            for qx = idcsQPU(ex) : idcsQPU(ex + 1) - 1
                k = 1;
                for j = 1 : NURBS.Order(2) + 1
                    for i = 1 : NURBS.Order(1) + 1
                        N0(k) = Nx(qx, i, 1) * Ny(qy, j, 1);
                        N1(1, k) = Nx(qx, i, 2) * Ny(qy, j, 1);
                        N1(2, k) = Nx(qx, i, 1) * Ny(qy, j, 2);
                        k = k + 1;
                    end
                end
                [~, R1] = Rationalize(Weights(Mesh.El(e, :)), N0, N1);
                W = Wx(qx) * Wy(qy);
                
                % Gradient of mapping from parameter space to physical space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                
                J1 = abs(det(dxdxi));
                
                %Compute derivatives of basis functions w.r.t physical coordinates
                dRdx = dxdxi \ R1;
                
                % B matrix
                B = zeros(3, 2 * Mesh.NEN);
                
                %B(1, 1:2:2*Mesh.NEN) = dRdx(1, :);
                %B(2, 2:2:2*Mesh.NEN) = dRdx(2, :);
                %B(3, 1:2:2*Mesh.NEN) = dRdx(2, :);
                %B(3, 2:2:2*Mesh.NEN ) = dRdx(1, :);
                
                % way to sort dispacelement: [u1 u2.....u_Ndof v1,v2,...]
                B(1, 1 : Mesh.NEN) = dRdx(1, :);
                B(2, Mesh.NEN + 1 : end) = dRdx(2, :);
                B(3, 1 : Mesh.NEN) = dRdx(2, :);
                B(3, Mesh.NEN + 1 : end) = dRdx(1, :);
                
                % compute element B and WeJ at quadrature point
                Btemp = zeros(3, Mesh.NDof);
                Btemp(:, getDofNEl) = B;
                countGP = countGP + 1;
                Bstrain{countGP} = sparse(Btemp);
                WeJ(countGP) = J1*W;
            end
        end
    end
% else
%     idcsX = findAllSpans(NURBS.NCtrlPts(1), NURBS.Order(1), NURBS.KntVect{1}, Mesh.NElDir(1));
%     idcsY = findAllSpans(NURBS.NCtrlPts(2), NURBS.Order(2), NURBS.KntVect{2}, Mesh.NElDir(2));
%     
%     % calll quaddrature rule for computing enhanced strain
%     [Xg, W] = GaussRule(NGPs(1));
%     [Yg, W] = GaussRule(NGPs(2));
%     
%     % add three incompatible modes
%     imode = 3;
%     
%     for ey = 1 : Mesh.NElDir(2)
%         for ex = 1 : Mesh.NElDir(1)
%             e = sub2ind(Mesh.NElDir, ex, ey);
%             
%             sctrEas = [Mesh.El(e, :), Mesh.NDof / Mesh.Dof + imode*e - 2, Mesh.NDof / Mesh.Dof + imode*e - 1, Mesh.NDof / Mesh.Dof + imode*e];
%             getDofNEl(1 : Mesh.NEN + imode) = sctrEas;
%             getDofNEl(Mesh.NEN + imode + 1 : 2*Mesh.NEN + 2*imode) = sctrEas + Mesh.NDof / Mesh.Dof + Mesh.NEl*imode;
%             
%             % get physical coordinates of knot vector
%             coordPhys = [];
%             for kc = 1 : length(Q4) % 4 vertexes of Q4
%                 u_hat = Q4(kc, 1); % quadrature point
%                 v_hat = Q4(kc, 2); % quadrature point
%                 Xi = (u_hat + 1) * ((NURBS.KntVect{1}(idcsX(ex) + 1) - NURBS.KntVect{1}(idcsX(ex))) / 2) + NURBS.KntVect{1}(idcsX(ex));
%                 Eta = (v_hat + 1) * ((NURBS.KntVect{2}(idcsY(ey) + 1) - NURBS.KntVect{2}(idcsY(ey))) / 2) + NURBS.KntVect{2}(idcsY(ey));
%                 tmp = BsplineEval(NURBS.KntVect, NURBS.CtrlPts4D, {Xi, Eta});
%                 x_physis = (tmp(1 : 3, :) ./ tmp(4, :))';
%                 coordPhys = [coordPhys; x_physis]; % gcoord contans coordinates of control points on each element (or patch)
%                 plot(coordPhys(:, 1), coordPhys(:, 2), 'kx'); hold on;
%             end
%             %             plot(coordPhys(:, 1), coordPhys(:, 2), 'kx'); hold on;
%             
%             % compute displacement strain matrix
%             for qy = 1 : NGPs(2)
%                 v_hat=Yg(qy);
%                 for qx = 1 : NGPs(1)
%                     u_hat=Xg(qx);
%                     k = 1;
%                     for j = 1 : NURBS.Order(2) + 1
%                         for i = 1 : NURBS.Order(1) + 1
%                             N0(k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 1);
%                             N1(1, k) = Nx(ex, qx, i, 2) * Ny(ey, qy, j, 1);
%                             N1(2, k) = Nx(ex, qx, i, 1) * Ny(ey, qy, j, 2);
%                             k = k + 1;
%                         end
%                     end
%                     % evalute EAS strains for incompressible modes
%                     dNdxip = calcBFemQ4([0, 0]);
%                     %                     dNdxip = calcBFemQ4([u_hat, v_hat]);             % Q4 element shape functions
%                     Jp = dNdxip'*coordPhys(:, 1 : 2);                  % element Jacobian matrix
%                     
%                     dMdxi12 = [-2*u_hat*(1 - v_hat^2), -2*v_hat*(1 - u_hat^2)];
%                     dMdxi34 = [-2*(1 - v_hat^2), -2*(1 - u_hat^2)];
%                     dMdxi56 = [4*u_hat 4*v_hat];
%                     
%                     Ba12 = [dMdxi12(2), dMdxi12(1);
%                         dMdxi34(1), dMdxi34(2)] / Jp;
%                     Ba56 = [dMdxi56(1), dMdxi56(2)] / Jp;
%                     
%                     % B_alpha = [Ba12(1,1)    Ba12(2,1)    Ba56(1)   0          0          0;
%                     %             0            0            0         Ba12(1,2)  Ba12(2,2)  Ba56(2);
%                     %             Ba12(1,2)    Ba12(2,2)    0         Ba12(1,1)  Ba12(2,1)   0];
%                     
%                     % Gradient of NURBS basis function
%                     [R0, R1] = Rationalize(Weights(Mesh.El(e, :)), N0, N1);
%                     J2 = Jx(ex) * Jy(ey);
%                     W = Wx(qx) * Wy(qy);
%                     
%                     % Gradient of mapping from parameter space to physical space
%                     dxdxi = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
%                     
%                     X = R0 * CtrlPts(Mesh.El(e, :), 1 : 2);
%                     
%                     plot(X(1), X(2), 'bh'); hold on;
%                     
%                     J1 = abs(det(dxdxi));
%                     
%                     %Compute derivatives of basis functions w.r.t physical coordinates
%                     dRdx = dxdxi \ R1;
%                     
%                     % B matrix
%                     B = zeros(3, 2 * Mesh.NEN + 2 * imode);
%                     B(1, 1 : Mesh.NEN + imode) = [dRdx(1, :), Ba12(1, 1), Ba12(2, 1), Ba56(1)];
%                     B(2, Mesh.NEN + imode + 1 : end) = [dRdx(2, :), Ba12(1, 2), Ba12(2, 2), Ba56(2)];
%                     B(3, 1 : Mesh.NEN + imode) = [dRdx(2, :), Ba12(1, 2), Ba12(2, 2), 0];
%                     B(3, Mesh.NEN + imode + 1 : end) = [dRdx(1, :), Ba12(1, 1), Ba12(2, 1), 0];
%                     
%                     % compute element B and WeJ at quadrature point
%                     Btemp = zeros(3, Mesh.NDof + 2*Mesh.NEl*imode);
%                     Btemp(:, getDofNEl) = B;
%                     countGP = countGP + 1;
%                     Bstrain{countGP} = sparse(Btemp);
%                     WeJ(countGP) = J1 * J2 * W;
%                 end
%             end
%         end
%     end
% end

%end
