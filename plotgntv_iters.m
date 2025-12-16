close all;clear; clc;
load('P_tv_results','P_all_bef')
% load('bsc_params_JCCP_GNTV.mat','params_all_reg');

fnames={'alo_smerino_iters/';
    'alo_smerino_iters_aniso/';
    'alo_smerino_iters_iso/';
'alo_smerino_iters_soloan/';
'alo_smerino_iters_fminsearch/';
'alo_smerino_iters_fminsearch_soloan/'};

fi=1;

lstfiles = ls([fnames{fi},'*.mat']);
lstfiles = lstfiles(1:6,:);

ncol=3;
nrow=ceil((size(lstfiles,1)+1)/ncol);

allMsize = [];
allMsize_lastscan = [];

for i=1:size(lstfiles,1)
    lblfile = strtrim(lstfiles(i,:));
    lambda = lblfile(29:end-4);
    load([fnames{fi},lblfile],'params_all_reg');
    cpFields = fieldnames(params_all_reg);
    M = [];
    
    for c = 1:numel(cpFields)
        scanFields = fieldnames(params_all_reg.(cpFields{c}));
        for s = 1:numel(scanFields)
            data = params_all_reg.(cpFields{c}).(scanFields{s});  % e.g., [52×4×7]
            
            % Reshape each scan to [N × 7]
            reshaped = reshape(data, [], 7);  % flatten first two dims
            
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
    plot(M(:,1)*1e6,M(:,2),'+b','DisplayName','N');hold on;
    plot(M(:,4)*1e6,M(:,5),'+r','DisplayName','C');
    xlim([0 12]);ylim([0 1]);
    xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
    title(['\lambda=',lambda]);
    xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
    xline(9,'k--','HandleVisibility','off','LineWidth',2);
    yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
    yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
    plot(mean(M(:,1)*1e6),mean(M(:,2)),...
        Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
    plot(mean(M(:,4)*1e6),mean(M(:,5)),...
        Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
    legend; axis square;theme light;

    figure(2);
    subplot(nrow,ncol,i+1);
    imagesc(data(:,:,1).');
    title(['\lambda=',lambda]); sgtitle([cpFields{c}, ' ', scanFields{s}]);
    colorbar; clim([0 12]*1e-6);

end

lastM = min(allMsize);
lastM_lastscan = min(allMsize_lastscan);

figure(1);
sgtitle('lsqnonlin')
% sgtitle('fminsearch')
subplot(nrow,ncol,1);
plot(P_all_bef(1:lastM,1)*1e6,P_all_bef(1:lastM,2),'+b','DisplayName','N');hold on;
plot(P_all_bef(1:lastM,4)*1e6,P_all_bef(1:lastM,5),'+r','DisplayName','C');
xlim([0 12]);ylim([0 1]);
xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
title('before TV');
xline(5.49,'k--','HandleVisibility','off','LineWidth',2);
xline(8.73,'k--','HandleVisibility','off','LineWidth',2);
yline(0.13,'k--','HandleVisibility','off','LineWidth',2);
yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
plot(mean(P_all_bef(1:lastM,1)*1e6),mean(P_all_bef(1:lastM,2)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
plot(mean(P_all_bef(1:lastM,4)*1e6),mean(P_all_bef(1:lastM,5)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
legend; axis square; theme light;

map0 = nan(size(data(:,:,1)));
map0_ind = find(~isnan(data(:,:,1)));
map0(map0_ind) = P_all_bef(1:lastM_lastscan,1);

figure(2);
sgtitle('lsqnonlin')
% sgtitle('fminsearch')
subplot(nrow,ncol,1);
imagesc(map0.'); title('before TV'); colorbar; clim([0 12]*1e-6);