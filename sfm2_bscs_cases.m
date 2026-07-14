clear; clc;
addpath('utils/');

% Each uncommented row is processed.
% Comment out any complete row that you do not want to run.
% Columns: case label, input MAT file, BSC variable, output MAT file
cases = {
    % 'JC',   'data/bscdataJCCP.mat',   'bscdataJCCP',   'data/sfm2_bsc_params_JCCP.mat';
    % 'LMTK', 'data/bscdataLMTKCP.mat', 'bscdataLMTKCP', 'data/sfm2_bsc_params_LMTKCP.mat';
    '4T1',  'data/bscdata4T1CP.mat',  'bscdata4T1CP',  'data/sfm2_bsc_params_4T1CP.mat';
};

totalTimer = tic;

for caseIdx = 1:size(cases, 1)
    caseTimer = tic;

    caseName   = cases{caseIdx, 1};
    inputFile  = cases{caseIdx, 2};
    dataVar    = cases{caseIdx, 3};
    outputFile = cases{caseIdx, 4};

    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('CASE %d/%d: %s\n', caseIdx, size(cases, 1), caseName);
    fprintf('Input : %s\n', inputFile);
    fprintf('Output: %s\n', outputFile);
    fprintf('%s\n', repmat('=', 1, 70));

    loadedData = load(inputFile);
    bscdata = loadedData.(dataVar);
    f = loadedData.f;

    cpNames = fieldnames(bscdata);
    params_all = struct();

    for c = 1:numel(cpNames)
        cpName = cpNames{c};
        scanNames = fieldnames(bscdata.(cpName));

        fprintf('\n  Sample %d/%d: %s (%d scans)\n', c, numel(cpNames), cpName, numel(scanNames));

        for s = 1:numel(scanNames)
            scanTimer = tic;
            scanName = scanNames{s};

            bscblock = bscdata.(cpName).(scanName);

            if isempty(bscblock)
                fprintf('    [%2d/%2d] %-20s SKIPPED — empty block\n', s, numel(scanNames), scanName);
                continue;
            end

            validMask = any(bscblock ~= 0, 3);
            [rows, cols] = find(validMask);
            nValid = numel(rows);

            fprintf('    [%2d/%2d] %-20s %7d pixels ... ', s, numel(scanNames), scanName, nValid);

            if nValid == 0
                fprintf('SKIPPED — no valid pixels\n');
                continue;
            end

            params_block = nan(size(bscblock, 1), size(bscblock, 2), 11);
            params_valid = nan(nValid, 11);

            parfor k = 1:nValid
                x = rows(k);
                y = cols(k);

                bsc_vector = squeeze(bscblock(x, y, :));

                params_valid(k, :) = sfm2_inversion_BSC_SFM_Neldermead_sansLog_Fc(f, bsc_vector, 20);
            end

            for k = 1:nValid
                params_block(rows(k), cols(k), :) = params_valid(k, :);
            end

            params_all.(cpName).(scanName) = params_block;

            fprintf('DONE in %.1f s\n', toc(scanTimer));
        end
    end

    save(outputFile, 'params_all', 'f', '-v7.3');

    fprintf('\nCompleted %s in %.1f min\n', caseName, toc(caseTimer) / 60);
end

fprintf('\nAll selected cases completed in %.1f min.\n', toc(totalTimer) / 60);