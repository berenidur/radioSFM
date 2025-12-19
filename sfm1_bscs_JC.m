clear; clc;

load('bscdataJCCP.mat'); % loads bscdataJCCP and f
f=f*1e6;
cpNames = fieldnames(bscdataJCCP);
params_all = struct();

for c = 1:numel(cpNames)
    cpName = cpNames{c};
    fprintf('\n=== Processing %s ===\n', cpName);
    scanNames = fieldnames(bscdataJCCP.(cpName));

    for s = 1:numel(scanNames)
        scanName = scanNames{s};
        fprintf('   -> %s\n', scanName);

        % Access the 3D BSC matrix: [x, y, frequency]
        bscblock = bscdataJCCP.(cpName).(scanName);

        % Skip if empty
        if isempty(bscblock)
            continue;
        end

        % Mask of pixels that are not entirely zero along 3rd axis
        % (true if any frequency value is nonzero)
        validMask = any(bscblock ~= 0, 3);

        % Get coordinates of valid pixels
        [rows, cols] = find(validMask);
        nValid = numel(rows);
        % fprintf('      Nonzero pixels: %d\n', nValid);

        params_block = nan(size(bscblock,1),size(bscblock,2),11);

        for k = 1:nValid
            x = rows(k);
            y = cols(k);

            bsc_vector = squeeze(bscblock(x, y, :));  % (128×1)

            params = sfm1_inversion_BSC_SFM_Neldermead_sansLog_Fc(f, bsc_vector, 20);
            
            params_block(x,y,:) = params;

        end

        params_all.(cpName).(scanName) = params_block;
    end
end

save('sfm1_bsc_params_JCCP.mat','params_all');