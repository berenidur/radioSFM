clear; clc;

%% === Load data ===
load('bscdataJCCP.mat'); % contains struct bscdataJCCP.(cpName).(scanName)
load('bsc_f.mat');       % vector f (frequencies)
load('bsc_params_JCCP.mat','params_all'); % your existing parameter maps
addpath('utils/')

cpNames = fieldnames(bscdataJCCP);
params_all_reg = struct();

for c = 1:numel(cpNames)
    cpName = cpNames{c};
    scanNames = fieldnames(bscdataJCCP.(cpName));

    for s = 1:numel(scanNames)
        scanName = scanNames{s};
        fprintf('\n=== Regularizing %s / %s ===\n', cpName, scanName);

        bscblock = bscdataJCCP.(cpName).(scanName);
        if isempty(bscblock)
            continue;
        end

        params_block = params_all.(cpName).(scanName);
        params_block = params_block(:,:,1:7);

        params_block(:,:,[1,4]) = 1e5 * params_block(:,:,[1,4]);

        %% Build mask of valid pixels (same as your original code)
        validMask = any(bscblock ~= 0, 3);
        if ~any(validMask(:))
            fprintf('No valid pixels, skipping.\n');
            continue;
        end

        %% === Configure GNTV parameters ===
        cfg.channels = 1:7;                 % parameters to regularize
        cfg.lambda   = 0.03 * ones(1,length(cfg.channels));    % TV weights
        cfg.rho      = 1;                   % ADMM penalty
        cfg.maxIter  = 10;                  % outer ADMM iterations
        cfg.tol      = 1e-4;
        cfg.lsq_opts = optimoptions('lsqnonlin', ...
            'Display','off', ...
            'Algorithm','levenberg-marquardt', ...
            'MaxIterations',100, ...
            'StepTolerance',1e-6, ...
            'FunctionTolerance',1e-6);
        cfg.display = true;
        cfg.validMask = validMask;          % << add mask to config

        %% === Run GNTV regularization ===
        params_block_reg = gntv_lsqnonlin(params_block, bscblock, f, cfg);

        %% === Reinstate NaNs outside valid region ===
        for ch = cfg.channels
            tmp = params_block_reg(:,:,ch);
            tmp(~validMask) = NaN;
            params_block_reg(:,:,ch) = tmp;
        end

        %% Save
        params_block_reg(:,:, [1,4]) = params_block_reg(:,:, [1,4]) / 1e5;
        params_all_reg.(cpName).(scanName) = params_block_reg;
        save('bsc_params_JCCP_GNTV.mat','params_all_reg','-v7.3');
    end
end

save('bsc_params_JCCP_GNTV.mat','params_all_reg','-v7.3');
fprintf('\nSaved: bsc_params_JCCP_GNTV.mat\n');
