function S = SPThickWalledCylindExactStress(Ri, Ro, pi, po, nu)
% S = SPThickWalledCylindExactStress(Ri, Ro, pi, po, nu)
% Evaluate analytical solution of stresses for "Thick-Walled Cylinder"
% problem
% ------------------------------------------------------------------
% Input:
%       Ri: inner radius
%       Ro: outer radius
%       pi: inner pressure
%       po: outer pressure
%       nu: Poisson's ratio
% ------------------------------------------------------------------
% Output:
%       S: a structure contains 4 stress functions,
%           S.xx(x, y)
%           S.yy(x, y)
%           S.xy(x, y)
%           S.zz(x, y)
% ------------------------------------------------------------------

r = @(x, y) sqrt(x .^ 2 + y .^ 2);
theta =  @(x, y) atan(y ./ x);
C1t =  @(x, y) cos(theta(x, y));
S1t =  @(x, y) sin(theta(x, y));

term1 = (pi * Ri ^ 2 - po * Ro ^ 2);
term2 = (Ro ^ 2 - Ri ^ 2);
term3 = Ri ^ 2 * Ro ^ 2 * (po - pi);
term4 =  @(x, y) r(x, y) .^ 2 .* (Ro ^ 2 - Ri ^ 2);

srr =  @(x, y) term1 / term2 + term3 ./ term4(x, y); % \sigma_{rr}
stt =  @(x, y) term1 / term2 - term3 ./ term4(x, y); % \sigma_{\theta \theta}
srt = 0;
szz =  @(x, y) nu * (srr(x, y) + stt(x, y)); % \sigma_{zz}

S.xx = @(x, y)srr(x, y) .* C1t(x, y) .^ 2 - 2 .* S1t(x, y) .* C1t(x, y) .* srt - stt(x, y) .* C1t(x, y) .^ 2 + stt(x, y); % \sigma_{xx}
S.yy = @(x, y)-srr(x, y) .* C1t(x, y) .^ 2 + 2 .* S1t(x, y) .* C1t(x, y) .* srt + C1t(x, y) .^ 2 .* stt(x, y) + srr(x, y); % \sigma_{yy}
S.xy = @(x, y) S1t(x, y) .* C1t(x, y) .* srr(x, y) - S1t(x, y) .* C1t(x, y) .* stt(x, y) + 2 .* C1t(x, y) .^ 2 .* srt - srt; % % \sigma_{xy}
S.zz = @(x, y) szz(x, y);
end



