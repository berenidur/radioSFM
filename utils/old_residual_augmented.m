function r = residual_augmented(x, f, data_vec, rho, u_tilde)
% Residuals = [data misfit; sqrt(rho)*(x - u_tilde)]
kk = 2*pi*f/1540;
bsc_model = gntv_sfm2_F_compute_bsc_from_x([x(:)], kk);
data_vec(isnan(data_vec)) = 0;

r_data = data_vec(:) - bsc_model(:);
r_reg  = sqrt(rho) * (x(:) - u_tilde(:));
r = [r_data; r_reg];
end
