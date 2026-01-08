close all;clear; clc;
% cp='JC';    hist = [5.49, 0.13, 8.73, 0.52]; % a_n, phi_n, a_c, phi_c (a in um)
cp='4T1';   hist = [5.35, 0.14, 8.12, 0.49];
% cp='LMTK';  hist = [4.66, 0.15, 7.25, 0.58];

% sfm = 1;
sfm = 2;

load(['data/sfm',num2str(sfm),'_bsc_params_',cp,'CP.mat'],'params_all')

fnames={'alo_smerino_iters/';               % 1
    'alo_smerino_iters_aniso/';             % 2
    'alo_smerino_iters_iso/';               % 3
'alo_smerino_iters_soloan/';                % 4
'alo_smerino_iters_fminsearch/';            % 5
'alo_smerino_iters_fminsearch_soloan/'};    % 6

fi=3;
% sfm2 iso sale bien

lstfiles = ls([fnames{fi},['sfm',num2str(sfm),'*',cp,'*.mat']]);
disp(lstfiles)
% keyboard;

ncol=3;
nrow=ceil((size(lstfiles,1)+1)/ncol);

allMsize = [];
allMsize_lastscan = [];

for i=1:size(lstfiles,1)
    lblfile = strtrim(lstfiles(i,:)); %disp(lblfile);    
    lambda = extractBetween(lblfile, 'lambda_', '.mat');
    load([fnames{fi},lblfile],'params_all_reg');
    cpFields = fieldnames(params_all_reg);
    M = [];
    
    for c = 1:numel(cpFields)
        scanFields = fieldnames(params_all_reg.(cpFields{c}));
        for s = 1:numel(scanFields)
            data = params_all_reg.(cpFields{c}).(scanFields{s});  % e.g., [52×4×7]
            
            % Reshape each scan to [N × 7]
            reshaped = reshape(data, [], -1+sfm*4);  % flatten first two dims
            
            % Concatenate along the first dimension
            M = [M; reshaped];
        end
    end
    M(any(isnan(M), 2), :) = [];
    reshaped(any(isnan(reshaped), 2), :) = [];
    
    allMsize = [allMsize, size(M,1)];
    allMsize_lastscan = [allMsize_lastscan, size(reshaped,1)];
    
    
    figure(1);
    subplot(nrow,ncol,i+1);
    if sfm==2
        plot(M(:,1)*1e6,M(:,2),'+b','DisplayName','N');hold on;
        plot(M(:,4)*1e6,M(:,5),'+r','DisplayName','C');
    elseif sfm==1
        plot(M(:,1)*1e6,M(:,2),'+b','DisplayName','N,C');hold on;
    end
    xlim([0 12]);ylim([0 1]);
    xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
    title(['\lambda=',lambda]);
    xline(hist(1),'k--','HandleVisibility','off','LineWidth',2);
    xline(hist(3),'k--','HandleVisibility','off','LineWidth',2);
    yline(hist(2),'k--','HandleVisibility','off','LineWidth',2);
    yline(hist(4),'k--','HandleVisibility','off','LineWidth',2);
    plot(mean(M(:,1)*1e6),mean(M(:,2)),...
        Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
    if sfm==2
        plot(mean(M(:,4)*1e6),mean(M(:,5)),...
            Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
    end
    legend; axis square;theme light;

    figure(2);
    subplot(nrow,ncol,i+1);
    imagesc(data(:,:,1).');
    title(['\lambda=',lambda]); sgtitle([cpFields{c}, ' ', scanFields{s}]);
    colorbar; clim([0 12]*1e-6);

end

lastM = min(allMsize);
lastM_lastscan = min(allMsize_lastscan);


fields = fieldnames(params_all.CP1);
C = cell(numel(fields),1);

for k = 1:numel(fields)
    A = params_all.CP1.(fields{k});
    X = reshape(A, [], size(A,3));
    C{k} = X(~isnan(X(:,1)), :);
end

before_reg = vertcat(C{:});

figure(1);
sgtitle([cp,' - lsqnonlin'])
% sgtitle('fminsearch')
subplot(nrow,ncol,1);
plot(before_reg(1:lastM,1)*1e6,before_reg(1:lastM,2),'+b','DisplayName','N');hold on;
plot(before_reg(1:lastM,4)*1e6,before_reg(1:lastM,5),'+r','DisplayName','C');
xlim([0 12]);ylim([0 1]);
xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
title('before TV');
xline(hist(1),'k--','HandleVisibility','off','LineWidth',2);
xline(hist(3),'k--','HandleVisibility','off','LineWidth',2);
yline(hist(2),'k--','HandleVisibility','off','LineWidth',2);
yline(hist(4),'k--','HandleVisibility','off','LineWidth',2);
plot(mean(before_reg(1:lastM,1)*1e6),mean(before_reg(1:lastM,2)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
plot(mean(before_reg(1:lastM,4)*1e6),mean(before_reg(1:lastM,5)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
legend; axis square; theme light;

map0 = nan(size(data(:,:,1)));
map0_ind = find(~isnan(data(:,:,1)));
map0(map0_ind) = before_reg(1:lastM_lastscan,1);

figure(2);
sgtitle('lsqnonlin')
% sgtitle('fminsearch')
subplot(nrow,ncol,1);
imagesc(map0.'); title('before TV'); colorbar; clim([0 12]*1e-6);