function r = residual_global_tv_prior_aniso(x, params_init, channels, ...
                                      bscblock, f, mask_valid, ...
                                      idx_valid, lambda, tv_eps)
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

    % === ANISOTROPIC TV (LS residual version) ===
    % TV_aniso = sum |u(i+1,j) - u(i,j)| + |u(i,j+1) - u(i,j)|
    %
    % LS solver requires r^2 = lambda * |delta|
    % -> r = sqrt(lambda) * |delta|^(1/2)
    
    r_tv_cell = cell(C,1);
    mask_valid_2d = mask_valid;
    
    for c = 1:C
        chan = channels(c);
        map_c = params_cur(:,:,chan);
    
        max_edges_per_c = 2*M*N;
        r_edges = zeros(max_edges_per_c,1);
        edge_count = 0;
    
        % Horizontal |u(i+1,j) - u(i,j)|
        for i = 1:M-1
            for j = 1:N
                if mask_valid_2d(i,j) && mask_valid_2d(i+1,j)
                    delta = map_c(i+1,j) - map_c(i,j);
                    edge_count = edge_count + 1;
                    r_edges(edge_count) = sqrt(lambda(c)) * sqrt(abs(delta));
                end
            end
        end
    
        % Vertical |u(i,j+1) - u(i,j)|
        for i = 1:M
            for j = 1:N-1
                if mask_valid_2d(i,j) && mask_valid_2d(i,j+1)
                    delta = map_c(i,j+1) - map_c(i,j);
                    edge_count = edge_count + 1;
                    r_edges(edge_count) = sqrt(lambda(c)) * sqrt(abs(delta));
                end
            end
        end
    
        r_tv_cell{c} = r_edges(1:edge_count);
    end
    
    r_tv = vertcat(r_tv_cell{:});


    % === Total residual ===
    r = [r_data; r_tv];
end
