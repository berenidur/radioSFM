function [F,gradF] = sfm2_F_myfun_SFMTroisParam_trr_dr(z,connu)

% Optimization parameters:
% z(1) = a_N
% z(2) = phi_N
% z(3) = gammaZ_N
% z(4) = a_C - a_N
% z(5) = phi_C
% z(6) = gammaZ_C
% z(7) = w

% ============================================================
% Objective value
% ============================================================

x = z;

% Convert back to physical SFM2 parameters
x(4) = z(1) + z(4);   % a_C = a_N + delta_a

% Call original objective function
F = sfm2_F_myfun_SFMTroisParam(x,connu);


% ============================================================
% Numerical gradient with respect to z
% ============================================================

gradF = zeros(size(z));

% Scales for each optimization parameter
scale = [ ...
    1e-6, ...   % a_N
    0.1,  ...   % phi_N
    0.01, ...   % gammaZ_N
    1e-6, ...   % a_C - a_N
    0.1,  ...   % phi_C
    0.01, ...   % gammaZ_C
    0.1];       % w


for ii = 1:length(z)

    % Relative finite-difference step
    h = 1e-5 * max(abs(z(ii)),scale(ii));

    zp = z;
    zm = z;

    zp(ii) = zp(ii) + h;
    zm(ii) = zm(ii) - h;


    % Convert to physical parameters
    xp = zp; xp(4) = zp(1) + zp(4);
    xm = zm; xm(4) = zm(1) + zm(4);


    % Evaluate original objective
    Fp = sfm2_F_myfun_SFMTroisParam(xp,connu);
    Fm = sfm2_F_myfun_SFMTroisParam(xm,connu);


    % Central finite difference
    gradF(ii) = (Fp-Fm)/(2*h);

end

end