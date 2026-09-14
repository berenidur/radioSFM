function [F,gradF] = sfm2_F_myfun_SFMTroisParam_trr(x,connu)

% ============================================================
% Objective value
% ============================================================
F = sfm2_F_myfun_SFMTroisParam(x,connu);

% ============================================================
% Numerical gradient for fmincon trust-region-reflective
% ============================================================
gradF = zeros(size(x));

% Scales for each SFM2 parameter
scale = [ ...
    1e-6, ...   % a_N
    0.1,  ...   % phi_N
    0.01, ...   % gammaZ_N
    1e-6, ...   % a_C
    0.1,  ...   % phi_C
    0.01, ...   % gammaZ_C
    0.1];       % w

for ii = 1:length(x)

    % Relative finite-difference step
    h = 1e-5 * max(abs(x(ii)),scale(ii));

    xp = x;
    xm = x;

    xp(ii) = xp(ii) + h;
    xm(ii) = xm(ii) - h;

    Fp = sfm2_F_myfun_SFMTroisParam(xp,connu);
    Fm = sfm2_F_myfun_SFMTroisParam(xm,connu);

    gradF(ii) = (Fp-Fm)/(2*h);

end

end