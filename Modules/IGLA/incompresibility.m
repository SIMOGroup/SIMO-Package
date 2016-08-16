function [AeqIn] = incompresibility(sdof,Bstrain,totalGP)

% incompresibility condition at every gauss point
AeqIn = sparse(totalGP,sdof);
for i = 1:totalGP
    Bx = Bstrain{i}(1,:); 
    By = Bstrain{i}(2,:);
    AeqIn(i,:) = Bx + By; 
   clear Bx By;
end
