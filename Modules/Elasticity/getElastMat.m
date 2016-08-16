function D = getElastMat(E, nu, lab)

% Evaluate elasticity matrix

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