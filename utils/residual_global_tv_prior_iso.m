function r = residual_global_tv_prior_iso(x, params_init, channels, ...
                                      bscblock, f, mask_valid, ...
                                      idx_valid, lambda)
% Build the full residual:
%   r = [ r_data ; r_tv ]
%
% r_data : data fidelity residuals (BSC model vs measured)
% r_tv   : TV residuals built from spatial finite differences.
%
% x : vector of length nValid*C (parameters on valid pixels + chosen channels)

    [M,N,P] = size(params_init);
    C = numel(channels);
    nValid = numel(idx_valid);
    F = size(bscblock,3);      % number of frequencies

    % --- Rebuild parameter maps from x ---
    u_valid = reshape(x, [nValid, C]);       % (nValid x C)
    u_full = zeros(M*N, C);
    u_full(idx_valid, :) = u_valid;

    % params_cur: MxNxP with updated channels
    params_cur = params_init;
    for c = 1:C
        chan = channels(c);
        tmp = params_cur(:,:,chan);
        tmp = tmp(:);
        tmp = reshape(tmp, M*N, 1);
        tmp = tmp(:);
        tmp(idx_valid) = u_valid(:,c);
        params_cur(:,:,chan) = reshape(tmp, M, N);
    end

    % === (1) Data residual ===
    % For each valid pixel (i,j):
    %   r_data_ij = data_vec - model_vec
    %
    r_data = zeros(nValid*F, 1);
    kk = 2*pi*f/1540;

    for k = 1:nValid
        idx = idx_valid(k);
        [i,j] = ind2sub([M,N], idx);

        data_vec = squeeze(bscblock(i,j,:));
        if all(~isfinite(data_vec)) || all(data_vec==0)
            r_data((k-1)*F+(1:F)) = 0;
            continue;
        end

        x_full = squeeze(params_cur(i,j,:))';   % 1xP
        bsc_model = gntv_sfm2_F_compute_bsc_from_x(x_full, kk);

        data_vec(isnan(data_vec)) = 0;
        r_data((k-1)*F+(1:F)) = data_vec(:) - bsc_model(:);
    end

    % === ISOTROPIC TV (LS residual version) ===
    % TV_iso(i,j) = sqrt( (dx)^2 + (dy)^2 )
    %
    % Least-squares form:
    %     r^2 = lambda * sqrt(dx^2 + dy^2)
    %  => r = sqrt(lambda) * (dx^2 + dy^2)^(1/4)
    
    r_tv_cell = cell(C,1);
    mask_valid_2d = mask_valid;
    
    for c = 1:C
        chan = channels(c);
        map_c = params_cur(:,:,chan);
    
        % At most (M*N) pixels contribute one isotropic term
        r_edges = zeros(M*N,1);
        edge_count = 0;
    
        for i = 1:M
            for j = 1:N
                if ~mask_valid_2d(i,j)
                    continue;
                end
    
                % Compute dx if inside bounds and neighbor valid
                dx_valid = (i < M) && mask_valid_2d(i+1,j);
                if dx_valid
                    dx = map_c(i+1,j) - map_c(i,j);
                else
                    dx = 0;
                end
    
                % Compute dy if inside bounds and neighbor valid
                dy_valid = (j < N) && mask_valid_2d(i,j+1);
                if dy_valid
                    dy = map_c(i,j+1) - map_c(i,j);
                else
                    dy = 0;
                end
    
                % Skip pixels that have no valid neighbors
                if ~(dx_valid || dy_valid)
                    continue;
                end
    
                % Isotropic gradient magnitude
                grad_mag = sqrt(dx^2 + dy^2);
    
                % LS residual
                edge_count = edge_count + 1;
                r_edges(edge_count) = sqrt(lambda(c)) * grad_mag^(1/4);
            end
        end
    
        r_tv_cell{c} = r_edges(1:edge_count);
    end
    
    r_tv = vertcat(r_tv_cell{:});


    % === Total residual ===
    r = [r_data; r_tv];
end
