function r = residual_augmented_full(x_sub, i, j, params_init, channels, f, data_vec, rho, u_tilde)
    % Build full 7-parameter vector
    x_full = params_init(i,j,:);    % original values
    x_full = x_full(:)';
    x_full(channels) = x_sub;       % replace only regularized ones
    
    % Forward model
    kk = 2*pi*f/1540;
    bsc_model = gntv_sfm2_F_compute_bsc_from_x(x_full, kk);

    % residuals
    data_vec(isnan(data_vec)) = 0;
    r_data = data_vec(:) - bsc_model(:);
    r_reg = sqrt(rho) * (x_sub(:) - u_tilde(:));

    r = [r_data; r_reg];
end