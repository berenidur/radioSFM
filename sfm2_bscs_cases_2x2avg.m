close all; clear; clc;
addpath('utils/');

% ============================================================
% 2x2 averaging parameter
% ============================================================
% Minimum number of valid pixels required inside each 2x2 box.
% Allowed values: 1, 2, 3, or 4.
%
% Example:
%   1 -> at least 1 of 4 pixels must be valid
%   2 -> at least 2 of 4 pixels must be valid
%   3 -> at least 3 of 4 pixels must be valid
%   4 -> all 4 pixels must be valid
%
minValidPixelsPerBox = 1;
avgbox = '2x2';

if ~ismember(minValidPixelsPerBox, 1:4)
    error('minValidPixelsPerBox must be an integer between 1 and 4.');
end


% Each uncommented row is processed.
% Comment out any complete row that you do not want to run.
% Columns:
%   case label
%   input MAT file
%   BSC variable
%   parameter output MAT file
names = {
    'JC',...
    'LMTK',...
    '4T1',...
    };

cp = strcat(names,'CP');

cases = [ ...
    names(:), ...
    strcat('data/bscdata', cp(:), '.mat'), ...
    strcat('bscdata', cp(:)), ...
    strcat(['data',avgbox,'/minValidPixelsPerBox', ...
        num2str(minValidPixelsPerBox),'/sfm2_bsc_params_'], ...
        cp(:), '.mat') ...
];


totalTimer = tic;

for caseIdx = 1:size(cases, 1)

    caseTimer = tic;

    caseName   = cases{caseIdx, 1};
    inputFile  = cases{caseIdx, 2};
    dataVar    = cases{caseIdx, 3};
    outputFile = cases{caseIdx, 4};


    % ========================================================
    % Output file for REDUCED / AVERAGED BSC data
    %
    % Example:
    %   data2x2/minValidPixelsPerBox1/bscdataJCCP.mat
    % ========================================================
    outputFolder = fileparts(outputFile);

    bscOutputFile = fullfile( ...
        outputFolder, ...
        [dataVar '.mat']);


    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('CASE %d/%d: %s\n', caseIdx, size(cases, 1), caseName);
    fprintf('Input       : %s\n', inputFile);
    fprintf('Params output: %s\n', outputFile);
    fprintf('BSC output   : %s\n', bscOutputFile);
    fprintf('%s\n', repmat('=', 1, 70));


    % ========================================================
    % Load original BSC data
    % ========================================================
    loadedData = load(inputFile);

    bscdata_original = loadedData.(dataVar);
    f = loadedData.f;


    cpNames = fieldnames(bscdata_original);

    % Structure containing inversion results
    params_all = struct();

    % ========================================================
    % Structure containing the new spatially averaged BSC data
    %
    % It has exactly the same hierarchy as the original:
    %
    %   bscdata_avg.(cpName).(scanName)
    %
    % but each BSC block has dimensions:
    %
    %   ceil(n/2) x ceil(m/2) x 128
    % ========================================================
    bscdata_avg = struct();


    for c = 1:numel(cpNames)

        cpName = cpNames{c};
        scanNames = fieldnames(bscdata_original.(cpName));

        fprintf('\n  Sample %d/%d: %s (%d scans)\n', ...
            c, numel(cpNames), cpName, numel(scanNames));


        for s = 1:numel(scanNames)

            scanTimer = tic;
            scanName = scanNames{s};


            % ========================================================
            % Original BSC block: n x m x 128
            % ========================================================
            bscblock_original = ...
                bscdata_original.(cpName).(scanName);


            if isempty(bscblock_original)

                fprintf('    [%2d/%2d] %-20s SKIPPED — empty block\n', ...
                    s, numel(scanNames), scanName);

                % Preserve empty scans in averaged BSC structure
                bscdata_avg.(cpName).(scanName) = bscblock_original;

                continue;
            end


            % ========================================================
            % ORIGINAL VALID MASK
            %
            % Pixel is valid when at least one BSC value is non-zero.
            %
            % Size:
            %   n x m
            % ========================================================
            validMask_original = any(bscblock_original ~= 0, 3);

            [n, m, nFreq] = size(bscblock_original);


            % ========================================================
            % NON-OVERLAPPING 2x2 REDUCTION
            %
            % Edge boxes are NOT discarded:
            %
            % Normal interior box: 2 x 2
            % Odd final row:       1 x 2
            % Odd final column:    2 x 1
            % Odd row + column:    1 x 1
            %
            % Output dimensions:
            %
            %   n2 = ceil(n/2)
            %   m2 = ceil(m/2)
            % ========================================================

            n2 = ceil(n / 2);
            m2 = ceil(m / 2);


            % Reduced BSC block.
            %
            % Invalid output boxes remain zero.
            bscblock = zeros( ...
                n2, m2, nFreq, ...
                'like', bscblock_original);


            % Reduced validity mask
            validMask = false(n2, m2);


            % Number of valid original pixels contributing
            % to each reduced box
            validCount = zeros(n2, m2, 'uint8');


            for i = 1:n2

                % Normally two rows.
                %
                % For odd n, the final box contains only
                % the final original row.
                r1 = 2*i - 1;
                r2 = min(2*i, n);
                rr = r1:r2;


                for j = 1:m2

                    % Normally two columns.
                    %
                    % For odd m, the final box contains only
                    % the final original column.
                    c1 = 2*j - 1;
                    c2 = min(2*j, m);
                    cc = c1:c2;


                    % ------------------------------------------------
                    % Validity mask for this spatial box
                    %
                    % Possible dimensions:
                    %
                    %   2 x 2
                    %   2 x 1
                    %   1 x 2
                    %   1 x 1
                    % ------------------------------------------------
                    boxMask = validMask_original(rr, cc);

                    nBoxValid = nnz(boxMask);

                    validCount(i, j) = nBoxValid;


                    % ------------------------------------------------
                    % Box is valid only if it contains at least
                    % minValidPixelsPerBox valid pixels
                    % ------------------------------------------------
                    if nBoxValid >= minValidPixelsPerBox

                        validMask(i, j) = true;


                        % Number of actual spatial pixels in
                        % this box: 1, 2, or 4
                        nPixelsInBox = numel(boxMask);


                        % Extract BSC data and reshape:
                        %
                        % rows x cols x 128
                        %
                        % becomes
                        %
                        % nPixelsInBox x 128
                        %
                        boxData = reshape( ...
                            bscblock_original(rr, cc, :), ...
                            nPixelsInBox, nFreq);


                        % Keep ONLY valid spatial pixels
                        validBoxData = ...
                            boxData(boxMask(:), :);


                        % Average only the valid pixels
                        meanBSC = mean(validBoxData, 1);


                        % Store averaged BSC
                        bscblock(i, j, :) = ...
                            reshape(meanBSC, 1, 1, nFreq);

                    end
                end
            end


            % ========================================================
            % SAVE REDUCED BSC BLOCK INTO NEW STRUCTURE
            %
            % This is done BEFORE checking nValid so that even a scan
            % containing zero valid reduced boxes is preserved.
            %
            % Invalid boxes have BSC = 0.
            % Valid boxes contain the mean BSC of the valid pixels
            % from their original 2x2 region.
            % ========================================================
            bscdata_avg.(cpName).(scanName) = bscblock;


            % ========================================================
            % At this point:
            %
            % size(bscblock)
            %       = n2 x m2 x 128
            %
            % size(validMask)
            %       = n2 x m2
            %
            % validCount(i,j)
            %       = number of valid original pixels in that box
            % ========================================================

            [rows, cols] = find(validMask);
            nValid = numel(rows);


            fprintf( ...
                '    [%2d/%2d] %-20s %7d valid 2x2 boxes (%dx%d -> %dx%d) ... ', ...
                s, numel(scanNames), scanName, ...
                nValid, n, m, n2, m2);


            if nValid == 0

                fprintf('SKIPPED inversion — no valid boxes\n');

                % Optional: preserve parameter matrix as all NaNs
                params_all.(cpName).(scanName) = ...
                    nan(n2, m2, 11);

                continue;
            end


            % ========================================================
            % SFM2 inversion on REDUCED BSC data
            % ========================================================

            % Result has spatial dimensions n2 x m2
            params_block = nan(n2, m2, 11);

            params_valid = nan(nValid, 11);


            for k = 1:nValid

                x = rows(k);
                y = cols(k);

                % Averaged 128-point BSC vector
                bsc_vector = squeeze( ...
                    bscblock(x, y, :));

                params_valid(k, :) = ...
                    sfm2_inversion_BSC_SFM_Neldermead_sansLog_Fc( ...
                    f, bsc_vector, 20);
            end


            % ========================================================
            % Put inversion results back into spatial matrix
            % ========================================================
            for k = 1:nValid

                params_block(rows(k), cols(k), :) = ...
                    params_valid(k, :);

            end


            params_all.(cpName).(scanName) = params_block;


            fprintf('DONE in %.1f s\n', toc(scanTimer));

        end
    end


    % ========================================================
    % Create output directory
    % ========================================================
    mkdir(outputFolder);


    % ========================================================
    % SAVE PARAMETER RESULTS
    % ========================================================
    save( ...
        outputFile, ...
        'params_all', ...
        'f', ...
        'minValidPixelsPerBox', ...
        'avgbox', ...
        '-v7.3');

    fprintf('\n  Saved parameter file as:\n');
    fprintf('    %s\n', outputFile);


    % ========================================================
    % SAVE REDUCED / AVERAGED BSC DATA
    %
    % Rename bscdata_avg to the ORIGINAL variable name.
    %
    % Example:
    %
    % input:
    %   bscdataJCCP
    %
    % saved reduced file:
    %   bscdataJCCP.mat
    %
    % containing variable:
    %   bscdataJCCP
    %
    % This makes the reduced dataset compatible with code that
    % expects the original variable naming convention.
    % ========================================================

    bscSaveStruct = struct();

    bscSaveStruct.(dataVar) = bscdata_avg;
    bscSaveStruct.f = f;
    bscSaveStruct.minValidPixelsPerBox = minValidPixelsPerBox;
    bscSaveStruct.avgbox = avgbox;

    save( ...
        bscOutputFile, ...
        '-struct', ...
        'bscSaveStruct', ...
        '-v7.3');


    fprintf('  Saved averaged BSC file as:\n');
    fprintf('    %s\n', bscOutputFile);


    fprintf('\nCompleted %s in %.1f min\n', ...
        caseName, toc(caseTimer) / 60);

end


fprintf('\nAll selected cases completed in %.1f min.\n', ...
    toc(totalTimer) / 60);
