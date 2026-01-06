clear; clc;

%% === Load data ===
load('data/bscdata4T1CP.mat'); % contains struct bscdata4T1CP.(cpName).(scanName) + freqs
load('data/sfm1_bsc_params_4T1CP.mat','params_all'); % existing parameter maps
addpath('utils/');
savefolder = 'alo_smerino_iters_aniso/';
mkdir(savefolder)

cpNames = fieldnames(bscdata4T1CP);
params_all_reg = struct();

mkdir(savefolder);

%% Iteration
for lambda_i_exp=-4:5
    lambda_i=10^lambda_i_exp;
    if lambda_i_exp==-4, lambda_i=0; end

    % for c = 1:numel(cpNames)
    for c = 1
        cpName    = cpNames{c};
        scanNames = fieldnames(bscdata4T1CP.(cpName));
        nScans    = numel(scanNames);

        % temporary storage for parfor results
        params_block_reg_all = cell(nScans,1);

        parfor s = 1:nScans
            scanName = scanNames{s};
            fprintf('\n=== Regularizing %s / %s ===\n', cpName, scanName);

            bscblock = bscdata4T1CP.(cpName).(scanName);
            if isempty(bscblock)
                % leave params_block_reg_all{s} empty
                continue;
            end

            params_block = params_all.(cpName).(scanName);
            params_block = params_block(:,:,1:3);

            params_block(:,:,1) = 1e5 * params_block(:,:,1);

            % Build mask of valid pixels
            validMask = any(bscblock ~= 0, 3);
            if ~any(validMask(:))
                continue;
            end

            %% === Configure GNTV parameters ===
            cfg = struct();
            cfg.channels = 1:3;                 % parameters to regularize
            cfg.lambda   = lambda_i * ones(1,length(cfg.channels));    % TV weights
            cfg.rho      = 1;                   % ADMM penalty
            cfg.maxIter  = 200;                  % outer ADMM iterations
            cfg.tol      = 1e-4;
            cfg.lsq_opts = optimoptions('lsqnonlin', ...
                'Display','off', ...
                'Algorithm','levenberg-marquardt', ...
                'MaxIterations',200, ...
                'StepTolerance',1e-6, ...
                'FunctionTolerance',1e-6);
            cfg.display   = true;
            cfg.validMask = validMask;
            cfg.cpName    = cpName;
            cfg.scanName  = scanName;

            %% === Run GNTV regularization ===
            cfg.tvform  = 'aniso';
            cfg.sfm     = 1;
            params_block_reg = gntv_lsqnonlin(params_block, bscblock, f, cfg);

            %% === Reinstate NaNs outside valid region ===
            for ch = cfg.channels
                tmp = params_block_reg(:,:,ch);
                tmp(~validMask) = NaN;
                params_block_reg(:,:,ch) = tmp;
            end

            %% Scale back before storing
            params_block_reg(:,:,1) = params_block_reg(:,:,1) / 1e5;

            % store in cell for this s
            params_block_reg_all{s} = params_block_reg;
        end  % end parfor

        % Now we're back in normal (serial) code, safe to build the struct
        for s = 1:nScans
            if ~isempty(params_block_reg_all{s})
                scanName = scanNames{s};
                params_all_reg.(cpName).(scanName) = params_block_reg_all{s};
            end
        end
    end

    % Save *outside* parfor, also safe
    save([savefolder,'sfm1_bsc_params_4T1CP_GNTV_lambda_',num2str(lambda_i),'.mat'], ...
         'params_all_reg','-v7.3');
    fprintf('\nSaved: sfm1_bsc_params_4T1CP_GNTV_lambda_%g.mat\n', lambda_i);
end
