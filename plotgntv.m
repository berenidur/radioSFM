close all; clear; clc;
load('P_tv_results')
load('bsc_params_JCCP_GNTV.mat','params_all_reg');

cpFields = fieldnames(params_all_reg);
M = [];

% for c = 1:numel(cpFields)
for c = 1:3
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

lastM=size(M,1);
%%
% figure;
% % subplot(1,3,1);
% plot(P_all_bef(:,1)*1e6,P_all_bef(:,2),'+b','DisplayName','N');hold on;
% plot(P_all_bef(:,4)*1e6,P_all_bef(:,5),'+r','DisplayName','C');
% xlim([0 12]);ylim([0 1]);theme light;
% xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
% title('before TV');
% xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
% xline(9,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
% plot(mean(P_all_bef(:,1)*1e6),mean(P_all_bef(:,2)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% plot(mean(P_all_bef(:,4)*1e6),mean(P_all_bef(:,5)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% legend; axis square;
% % 
% % subplot(1,3,2);
% % plot(P_all_aft(:,1)*1e6,P_all_aft(:,2),'+b','DisplayName','N');hold on;
% % plot(P_all_aft(:,4)*1e6,P_all_aft(:,5),'+r','DisplayName','C');
% % xlim([0 12]);ylim([0 1]);theme light;
% % xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
% % title('after TV');
% % xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
% % xline(9,'k--','HandleVisibility','off','LineWidth',2);
% % yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
% % yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
% % plot(mean(P_all_aft(:,1)*1e6),mean(P_all_aft(:,2)),...
% %     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% % plot(mean(P_all_aft(:,4)*1e6),mean(P_all_aft(:,5)),...
% %     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% % legend; axis square;
%%
% figure;
% subplot(1,3,1);
% plot(P_all_bef(1:lastM,1)*1e6,P_all_bef(1:lastM,2),'+b','DisplayName','N');hold on;
% plot(P_all_bef(1:lastM,4)*1e6,P_all_bef(1:lastM,5),'+r','DisplayName','C');
% xlim([0 12]);ylim([0 1]);theme light;
% xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
% title('before TV');
% xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
% xline(9,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
% plot(mean(P_all_bef(1:lastM,1)*1e6),mean(P_all_bef(1:lastM,2)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% plot(mean(P_all_bef(1:lastM,4)*1e6),mean(P_all_bef(1:lastM,5)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% legend; axis square;
% 
% subplot(1,3,2);
% plot(P_all_aft(1:lastM,1)*1e6,P_all_aft(1:lastM,2),'+b','DisplayName','N');hold on;
% plot(P_all_aft(1:lastM,4)*1e6,P_all_aft(1:lastM,5),'+r','DisplayName','C');
% xlim([0 12]);ylim([0 1]);theme light;
% xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
% title('after TV');
% xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
% xline(9,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
% plot(mean(P_all_aft(1:lastM,1)*1e6),mean(P_all_aft(1:lastM,2)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% plot(mean(P_all_aft(1:lastM,4)*1e6),mean(P_all_aft(1:lastM,5)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% legend; axis square;
% 
% 
% subplot(1,3,3);
% plot(M(:,1)*1e6,M(:,2),'+b','DisplayName','N');hold on;
% plot(M(:,4)*1e6,M(:,5),'+r','DisplayName','C');
% xlim([0 12]);ylim([0 1]);theme light;
% xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
% title('new reg');
% xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
% xline(9,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
% yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
% plot(mean(P_all_aft(:,1)*1e6),mean(P_all_aft(:,2)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% plot(mean(P_all_aft(:,4)*1e6),mean(P_all_aft(:,5)),...
%     Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
% legend; axis square;

%%
figure;
subplot(1,2,1);
plot(P_all_bef(1:lastM,1)*1e6,P_all_bef(1:lastM,2),'+b','DisplayName','N');hold on;
plot(P_all_bef(1:lastM,4)*1e6,P_all_bef(1:lastM,5),'+r','DisplayName','C');
xlim([0 12]);ylim([0 1]);theme light;
xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
title('before TV');
xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
xline(9,'k--','HandleVisibility','off','LineWidth',2);
yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
plot(mean(P_all_bef(1:lastM,1)*1e6),mean(P_all_bef(1:lastM,2)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
plot(mean(P_all_bef(1:lastM,4)*1e6),mean(P_all_bef(1:lastM,5)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
legend; axis square;

subplot(1,2,2);
plot(M(:,1)*1e6,M(:,2),'+b','DisplayName','N');hold on;
plot(M(:,4)*1e6,M(:,5),'+r','DisplayName','C');
xlim([0 12]);ylim([0 1]);theme light;
xlabel('Scatterer radius (um)');ylabel('Scatterer volume fraction');
title('after TV');
xline(5.5,'k--','HandleVisibility','off','LineWidth',2);
xline(9,'k--','HandleVisibility','off','LineWidth',2);
yline(0.52,'k--','HandleVisibility','off','LineWidth',2);
yline(0.1,'k--','HandleVisibility','off','LineWidth',2);
plot(mean(P_all_aft(:,1)*1e6),mean(P_all_aft(:,2)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
plot(mean(P_all_aft(:,4)*1e6),mean(P_all_aft(:,5)),...
    Marker=".",MarkerSize=20,HandleVisibility="off",Color='k');
legend; axis square;

%%
figure;
subplot(2,2,1);imagesc(data(:,:,1).');colorbar;clim([0 12]*1e-6)
subplot(2,2,2);imagesc(data(:,:,2).');colorbar;clim([0 1])
subplot(2,2,3);imagesc(data(:,:,4).');colorbar;clim([0 12]*1e-6)
subplot(2,2,4);imagesc(data(:,:,5).');colorbar;clim([0 1])