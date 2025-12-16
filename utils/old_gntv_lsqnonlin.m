function params_block_reg = gntv_lsqnonlin(params_init, bscblock, f, cfg)
% Regularized parameter estimation (GNTV-like, Merino et al. 2025)
% Mask-aware version

[M,N,P] = size(params_init);
channels = cfg.channels;
C = numel(channels);

% Mask handling
if isfield(cfg,'validMask')
    mask_valid = cfg.validMask;
else
    mask_valid = any(~isnan(params_init(:,:,1)),3);
end
idx_valid = find(mask_valid);
nValid = numel(idx_valid);

u = zeros(M*N, C);
for c = 1:C
    tmp = params_init(:,:,channels(c));
    tmp(isnan(tmp)) = 0;
    u(:,c) = tmp(:);
end

v = u;
v_prev = v;
w = zeros(size(u));

rho = cfg.rho;
lambda = cfg.lambda;

fprintf('GNTV on %d valid pixels (%.1f%% of map)\n', nValid, 100*nValid/(M*N));
% === Residual tracking ===
res_primal = zeros(cfg.maxIter,1);
res_dual   = zeros(cfg.maxIter,1);

for l = 1:cfg.maxIter
    fprintf('--- ADMM iteration %d ---\n', l);
    u_prev = u;
    %% === Residual computation ===
    noNaNsind = ~isnan(v);
    res_primal(l) = norm(u(noNaNsind) - v(noNaNsind));           % primal residual
    res_dual(l)   = norm(v(noNaNsind) - v_prev(noNaNsind));      % dual residual (requires storing v_prev)

    %% === u-update ===
    for k = 1:nValid
        idx = idx_valid(k);
        [i,j] = ind2sub([M,N], idx);
        data_vec = squeeze(bscblock(i,j,:));
        if all(data_vec==0), continue; end

        u_tilde = v(idx,:) - w(idx,:);
        x0 = u(idx,:);
        
        % === Build full 7-parameter vector inside the objective ===
        fun = @(x_sub) residual_augmented_full(x_sub, i, j, params_init, channels, f, data_vec, rho, u_tilde);
        
        [xopt,~,~,~] = lsqnonlin(fun, x0, [], [], cfg.lsq_opts);
        
        u(idx,:) = xopt;
    end

    %% === v-update (TV denoising within mask) ===
    v_prev = v;   % save previous v for dual residual
    uplusw = u + w;
    for c = 1:C
        tmp = reshape(uplusw(:,c), M, N);
        tmp(~mask_valid) = 0;
        v_map = tv_denoise_chambolle_masked(tmp, mask_valid, lambda(c)/rho, 150);
        v(:,c) = v_map(:);
    end

    %% === w-update ===
    w = w + (u - v);

    relchg = norm(u(:)-u_prev(:)) / max(1e-8,norm(u_prev(:)));
    fprintf(' relchg = %.3e\n', relchg);
    if relchg < cfg.tol, break; end
end

params_block_reg = params_init;
for c = 1:C
    tmp = reshape(u(:,c), M,N);
    params_block_reg(:,:,channels(c)) = tmp;
end

% fig = figure('Visible','off'); 
% semilogy(res_primal(1:l), 'LineWidth', 1.5); hold on;
% semilogy(res_dual(1:l), 'LineWidth', 1.5);
% xlabel('Iteration');
% ylabel('Residual (log scale)');
% legend('Primal residual','Dual residual');
% title('ADMM Residuals');
% grid on;
% 
% % Save figure
% fname_fig = sprintf('iters_solo_an/%s_%s_residual_plot_iter_lambda_%s', cfg.cpName, cfg.scanName, num2str(lambda));
% saveas(fig, [fname_fig,'.png']);
% savefig(fig, [fname_fig,'.fig']);
% 
% % Close figure
% close(fig);

end
