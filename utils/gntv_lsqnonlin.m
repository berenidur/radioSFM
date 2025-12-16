function params_block_reg = gntv_lsqnonlin(params_init, bscblock, f, cfg)
% Global TV-regularized lsqnonlin:
%
%   min_x  0.5 * ||A(x) - data||^2  +  sum_c lambda_c * TV( x_c )
%
% TV is computed from spatial finite differences of the parameter maps
% and added as extra residuals. There is NO separate TV iteration.
%
% params_init : M x N x P   (initial parameter maps)
% bscblock    : M x N x F   (BSC spectra per pixel)
% f           : 1 x F       (frequencies)
% cfg         : struct:
%   .channels   channels to regularize (e.g. 1:7)
%   .lambda     TV weight per channel (1 x C)
%   .lsq_opts   options for lsqnonlin
%   .validMask  logical MxN of pixels to include
%   .tv_eps     (optional) small epsilon for TV smoothing
%
% Output:
%   params_block_reg : M x N x P final parameter maps

    [M,N,P] = size(params_init);
    channels = cfg.channels(:)';        % row
    C = numel(channels);

    % Mask of valid pixels
    if isfield(cfg,'validMask') && ~isempty(cfg.validMask)
        mask_valid = cfg.validMask;
    else
        mask_valid = ~isnan(params_init(:,:,channels(1)));
    end
    idx_valid = find(mask_valid);
    nValid = numel(idx_valid);

    if nValid == 0
        warning('No valid pixels, returning params_init unchanged.');
        params_block_reg = params_init;
        return;
    end

    lambda = cfg.lambda(:)';   % 1 x C
    if numel(lambda) ~= C
        error('cfg.lambda must have one entry per channel in cfg.channels.');
    end

    if isfield(cfg,'display') && cfg.display
        fprintf('Global TV-lsqnonlin on %d valid pixels (%.1f%% of map)\n', ...
            nValid, 100*nValid/(M*N));
    end

    % --- Build initial optimization vector x0 ---
    % We only optimize parameters on valid pixels and chosen channels.
    % x0 is of size nValid*C (column vector).
    u0_valid = zeros(nValid, C);
    for c = 1:C
        tmp = params_init(:,:,channels(c));
        tmp = tmp(:);
        u0_valid(:,c) = tmp(idx_valid);
    end
    x0 = u0_valid(:);  % vectorized

    % --- Create handle to residual function for lsqnonlin ---
    switch lower(cfg.tvform)
        case 'iso'
            fun = @(x) residual_global_tv_prior_iso( ...
                        x, params_init, channels, ...
                        bscblock, f, mask_valid, ...
                        idx_valid, lambda);    
        case 'aniso'
            fun = @(x) residual_global_tv_prior_aniso( ...
                        x, params_init, channels, ...
                        bscblock, f, mask_valid, ...
                        idx_valid, lambda);
        otherwise
            error('Unknown cfg.tvform value: %s (expected ''iso'' or ''aniso'')', cfg.tvform);
    end

    % --- Call lsqnonlin (Levenberg–Marquardt ~ Gauss–Newton) ---
    x_opt = lsqnonlin(fun, x0, [], [], cfg.lsq_opts);

    % --- Put result back into 3D parameter block ---
    u_opt_valid = reshape(x_opt, [nValid, C]);

    params_block_reg = params_init;
    for c = 1:C
        chan = channels(c);
        tmp = params_block_reg(:,:,chan);
        tmp = tmp(:);
        tmp(idx_valid) = u_opt_valid(:,c);
        params_block_reg(:,:,chan) = reshape(tmp, M, N);
    end

    % Outside valid mask, keep original or set NaN as you do later
end
