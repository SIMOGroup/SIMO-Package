function KVals = calcLocalStiffnessMatrices2DPFEM(Mesh, NURBS, E, nu, lab)
% function KVals = calcLocalStiffnessMatrices2D(Mesh, NURBS, E, nu, lab)
% -------------------------------------------------------------------
% Calculate two dimensional local stiffness matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       E: Young's modulus
%       nu: Poisson's ratio
%       lab: a string to determine stress state
%             ex. 'PlaneStress', 'PlaneStrain'
% -------------------------------------------------------------------
% Output:
%       KVals: matrix stores local stiffness matrices
%           (size(KVals) = [(NEN * 2) ^ 2, NEL])
% -------------------------------------------------------------------

D = getElastMat(E, nu, lab);

KVals = zeros((Mesh.NEN * Mesh.Dof) ^ 2, Mesh.NEl);

Weights = reshape(NURBS.Weights, 1, []);
CtrlPts = reshape(NURBS.CtrlPts3D, 3, [])';

idcsX = findAllSpans(NURBS.NCtrlPts(1), NURBS.Order(1), NURBS.KntVect{1}, Mesh.NElDir(1));
idcsY = findAllSpans(NURBS.NCtrlPts(2), NURBS.Order(2), NURBS.KntVect{2}, Mesh.NElDir(2));

for ey = 1 : Mesh.NElDir(2)
    for ex = 1 : Mesh.NElDir(1)
        e = sub2ind(Mesh.NElDir, ex, ey);
        Ke = zeros(Mesh.NEN * 2);
        % evaluate parametric coordinates of corner points
        u = [NURBS.KntVect{1}(idcsX(ex)), NURBS.KntVect{1}(idcsX(ex) + 1)];
        v = [NURBS.KntVect{2}(idcsY(ey)), NURBS.KntVect{2}(idcsY(ey) + 1)];
        
        uc = 0.5*sum(u);
        vc = 0.5*sum(v);
        
        [qPts, qPtsb, W, Wb] = getQuadData(3);
        nQuad = numel(W);
        nbQuad = numel(Wb);
        
        cell = zeros(5, 2);
        cell(1, :) = [u(1), v(1)];
        cell(2, :) = [u(2), v(1)];
        cell(3, :) = [u(2), v(2)];
        cell(4, :) = [u(1), v(2)];
        cell(5, :) = [uc, vc];
        
        cellConn = [1, 2, 5; 2, 3, 5; 3, 4, 5; 4, 1, 5]';
        
        for sc = 1 : 4 % loop over subcells
            term1 = zeros(2, Mesh.NEN); % \sum_{b=1}^{3} \sum_{j=1}^{2} w_j^b N_I n^b, where the 1st and 2nd row are multiplied by n_x, n_y respectively
            term2 = zeros(2, Mesh.NEN); % \sum_{b=1}^{3} \sum_{j=1}^{2} w_j^b N_I xy n_x^b
            term3 = zeros(2, Mesh.NEN); % \sum_{b=1}^{3} \sum_{j=1}^{2} w_j^b N_I xy n_y^b
            term4 = zeros(1, Mesh.NEN);
            for bd = 1 : 3 % loop over triangle edges
                % Loop through quadrature points.
                if bd == 1
                    coords = cell(cellConn(1 : 2, sc), :);
                elseif bd == 2
                    coords = cell(cellConn(2 : 3, sc), :);
                elseif bd == 3
                    coords = cell(cellConn(3 : -2 : 1, sc), :);
                end
                l = norm(coords(2, :) - coords(1, :));
                nx = (coords(2, 2) - coords(1, 2)) / l;
                ny = -(coords(2, 1) - coords(1, 1)) / l;
                normal = [nx; ny];
                dNdxib = [-1, 1] / 2;
                J2b = norm(dNdxib * coords);
                xi = qPtsb;
                Nb = [1 - xi, 1 + xi] / 2;
                uvb = Nb * coords;
                for nbq = 1 : nbQuad
                    N0xb = BasisFuns(idcsX(ex), uvb(nbq, 1), NURBS.Order(1), NURBS.KntVect{1});
                    N0yb = BasisFuns(idcsY(ey), uvb(nbq, 2), NURBS.Order(2), NURBS.KntVect{2});
                    
                    N0b = calcOuterProduct(N0xb, N0yb); % [1, (p + 1) * (q + 1)]
                    
                    R0b = RationalizeMod(Weights(Mesh.El(e, :)'), N0b);
                    
                    tmp1 = Wb(nbq) * R0b * J2b; % w_j ^ b * N_I
                    tmp2 = bsxfun(@times, tmp1, normal); % w_j ^ b * N_I * normal
                    term1 = term1 + tmp2;
                    
                    term2 = term2 + bsxfun(@times, tmp2(1, :), uvb(nbq, :)'); % w_j ^ b * N_I * n_x * x_j and w_j ^ b * N_I * n_x * y_j
                    term3 = term3 + bsxfun(@times, tmp2(2, :), uvb(nbq, :)'); % w_j ^ b * N_I * n_y * x_j and w_j ^ b * N_I * n_y * y_j
                end
            end
            % loop over interior gaussian points
            xi = qPts(:, 1); eta = qPts(:, 2);
            N = [1 - xi - eta, xi, eta]; % T3
            dNdxi = [-1, -1; 1, 0; 0, 1];
            uv = N * cell(cellConn(:, sc), :); % Gaussian coordinates in the parameter space
            J2 =  det(dNdxi' * cell(cellConn(:, sc), :)); % jacobian of mapping from parent space to parametric space
            WW = zeros(nQuad);
            JW = 0.5 * J2 * W; % jacobian of mapping
            WW(1, :) = JW';
            WW(2 : 3, :) = (bsxfun(@times, JW, uv))';
            for nq = 1 : nQuad
                N0x = BasisFuns(idcsX(ex), uv(nq, 1), NURBS.Order(1), NURBS.KntVect{1});
                N0y = BasisFuns(idcsY(ey), uv(nq, 2), NURBS.Order(2), NURBS.KntVect{2});
                
                N0 = calcOuterProduct(N0x, N0y);
                
                R0 = RationalizeMod(Weights(Mesh.El(e, :)'), N0);
                
                term4 = term4 + JW(nq) .* R0;
            end
            fu = zeros(3, Mesh.NEN);
            fv = zeros(3, Mesh.NEN);
            
            fu(1, :) = term1(1, :);
            fu(2, :) = term2(1, :) - term4;
            fu(3, :) = term2(2, :);
            dNdu = WW \ fu; clear fu
            
            fv(1, :) = term1(2, :);
            fv(2, :) = term3(1, :);
            fv(3, :) = term3(2, :) - term4;
            dNdv = WW \ fv; clear fv
            for nqp = 1 : nQuad
                R1(1, :) = dNdu(nqp, :);
                R1(2, :) = dNdv(nqp, :);
                % Gradient of mapping from parameter space to physical space
                dxdxi = R1 * CtrlPts(Mesh.El(e, :), 1 : 2);
                J1 = abs(det(dxdxi));
                
                %Compute derivatives of basis functions w.r.t physical coordinates
                dRdx = dxdxi \ R1;
                % B matrix
                %        _                                      _
                %        |  N_1,x  N_2,x  ...      0      0  ... |
                %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
                %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
                %        -                                      -
                B = zeros(3, 2 * Mesh.NEN);
                B(1, 1 : Mesh.NEN) = dRdx(1, :);
                B(2, Mesh.NEN + 1 : end) = dRdx(2, :);
                B(3, 1 : Mesh.NEN) = dRdx(2, :);
                B(3, Mesh.NEN + 1 : end) = dRdx(1, :);
                % compute element stiffness at quadrature point
                Ke = Ke + B' * D * B * JW(nqp) * J1;
            end
        end
        KVals(:, e) = Ke(:);
    end
end
end