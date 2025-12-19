% Estimation de trois parametres
% x(1)=rayon diffuseur                  a1
% x(2)=concentration                    phi1
% x(3)=contraste d'impedance gammaZ     gammaz1

function [xminbest_min,a,exitflag_min]=sfm1_F_InversionSFMTroisParam_NelderMead_Fc(connu,nb_xini)

% a_test_ini1=linspace(4e-6,10e-6,nb_xini);
% phi_test_ini1=linspace(0,0.7,nb_xini);
% gammaZ_test_ini1=linspace(0.005,0.055,nb_xini);
% 
% a_ini1=a_test_ini1(randperm(length(a_test_ini1)));  
% phi_ini1=phi_test_ini1(randperm(length(phi_test_ini1)));
% gammaZ_ini1=gammaZ_test_ini1(randperm(length(gammaZ_test_ini1)));

% load('sfm1_params_ini.mat','a_ini1','phi_ini1','gammaZ_ini1');
% vars = whos;
% for k = 1:length(vars)
%     name = vars(k).name;
%     val = eval(name);
%     disp([name ' = ' mat2str(val, 17) ';'])
% end
% keyboard;

a_ini1 = [5.5789473684210524e-06 4.631578947368421e-06 8.7368421052631584e-06 7.1578947368421058e-06 5.8947368421052634e-06 9.0526315789473686e-06 8.1052631578947381e-06 6.2105263157894736e-06 5.2631578947368422e-06 9.6842105263157906e-06 8.4210526315789482e-06 3.9999999999999998e-06 4.3157894736842108e-06 4.9473684210526312e-06 9.3684210526315788e-06 7.7894736842105279e-06 6.5263157894736846e-06 1.0000000000000001e-05 6.8421052631578957e-06 7.473684210526316e-06];
gammaZ_ini1 = [0.041842105263157896 0.039210526315789473 0.015526315789473683 0.055 0.031315789473684207 0.049736842105263163 0.0076315789473684215 0.0050000000000000001 0.020789473684210531 0.018157894736842106 0.044473684210526311 0.052368421052631578 0.02342105263157895 0.036578947368421058 0.033947368421052636 0.026052631578947369 0.028684210526315791 0.010263157894736842 0.012894736842105264 0.047105263157894733];
phi_ini1 = [0.40526315789473683 0.69999999999999996 0.036842105263157891 0.47894736842105262 0.51578947368421046 0.33157894736842103 0.44210526315789467 0.62631578947368416 0.25789473684210523 0.29473684210526313 0 0.14736842105263157 0.58947368421052626 0.36842105263157893 0.22105263157894733 0.18421052631578946 0.55263157894736847 0.073684210526315783 0.11052631578947367 0.66315789473684206];


options = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolX',1e-100);

for ii=1:length(a_ini1)
    % xinit=[1e-6 0.2 0.05];
    xinit=[
        a_ini1(ii), phi_ini1(ii), gammaZ_ini1(ii), ...
        % a_ini2(ii), phi_ini2(ii), gammaZ_ini2(ii), ...
        % w_ini(ii),
        ];
    [xmin,fvalmin,exitflag]=fminsearch(@(x) sfm1_F_myfun_SFMTroisParam(x,connu),xinit,options);
    % res=[xmin(1).*1e6 xmin(2) xmin(3)]
    if xmin(3)<0
       fvalmin=50;
    end
    xminbest(:,ii)=xmin;
    fvalminbest(ii)=fvalmin;
    exitflagbest(ii)=exitflag;
end


% [xminbest(1,:).*1e6; xminbest(2,:); xminbest(3,:)]
% fvalminbest

[a b]=min(fvalminbest);
xminbest_min=xminbest(:,b);
exitflag_min=exitflagbest(b);

%%%%%%%%%%%%%%
% fig=figure;
% plot(xminbest(1,:)*1e6,xminbest(2,:),'+');
% xlabel('Scatterer radius');ylabel('scatterer volume fraction');
% xlim([0 12]);ylim([0 1]);
% fig.Theme = "light"; 
% % keyboard;

end