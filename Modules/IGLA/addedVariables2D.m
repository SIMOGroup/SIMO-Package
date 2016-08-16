function [A2] = addedVariables2D(sdof, k1, k2, Bstrain, totalGP, stressState)       % Plane Strain case)

A2 = sparse(k2, sdof);
count = 1;
if ( strcmp(stressState,'PlaneStress') )       % Plane Strain case
    for i = 1 : totalGP
        Bx = Bstrain{i}(1, :); % see how to store this value in evaluated shape function
        By = Bstrain{i}(2, :);
        Bxy = Bstrain{i}(3, :);
        A2(count,:) = 2/sqrt(3)*Bx + 1/sqrt(3)*By;
        A2(count + 1, :) = By;
        A2(count + 2, :) = 1 / sqrt(3)*Bxy;
        count = count + 3;
        clear Bx By Bxy;
    end
else
    for i = 1 : totalGP
        Bx = Bstrain{i}(1, :); % see how to store this value in evaluated shape function
        By = Bstrain{i}(2, :);
        Bxy = Bstrain{i}(3, :);
        A2(count, :) = Bx - By;
        A2(count + 1, :) = Bxy;
        count = count + 2;
        clear Bx By Bxy;
    end
end
A2 = [A2, sparse(k2, k1), -sparse(1 : k2, 1 : k2, ones(k2, 1))];
% A2 = sparse([[A2], [zeros(k2,k1)], [-eye(k2)]]);