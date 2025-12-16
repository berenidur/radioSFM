% Estimation de sept parametres
% x(1)=rayon diffuseur                  a1          N
% x(2)=concentration                    phi1
% x(3)=contraste d'impedance gammaZ     gammaz1
% x(4)=rayon diffuseur                  a2          C
% x(5)=concentration                    phi2
% x(6)=contraste d'impedance gammaZ     gammaz2
% x(7)=n to c scattering ratio          w

function [xminbest_min,a,exitflag_min]=sfm2_F_InversionSFMTroisParam_NelderMead_Fc(connu,nb_xini)

% a_ini=(rand(1,nb_xini)*9+4).*1e-6;  étude CP/T inversion param pertinents
% phi_ini=(rand(1,nb_xini)*0.7);
% gammaZ_ini=(rand(1,nb_xini)*0.05+0.005);

% a_test_ini1=linspace(1e-6,7e-6,nb_xini);
% phi_test_ini1=linspace(0,0.3,nb_xini);
% gammaZ_test_ini1=linspace(0.005,0.055,nb_xini);
% 
% a_test_ini2=linspace(5e-6,9e-6,nb_xini);
% phi_test_ini2=linspace(0.45,0.75,nb_xini);
% gammaZ_test_ini2=linspace(0.005,0.055,nb_xini);
% 
% w_test_ini=linspace(0,1,nb_xini);
% 
% a_ini1=a_test_ini1(randperm(length(a_test_ini1)));  
% phi_ini1=phi_test_ini1(randperm(length(phi_test_ini1)));
% gammaZ_ini1=gammaZ_test_ini1(randperm(length(gammaZ_test_ini1)));
% 
% a_ini2=a_test_ini2(randperm(length(a_test_ini2)));  
% phi_ini2=phi_test_ini2(randperm(length(phi_test_ini2)));
% gammaZ_ini2=gammaZ_test_ini2(randperm(length(gammaZ_test_ini2)));
% 
% w_ini=w_test_ini(randperm(length(w_test_ini)));

% load('params_ini.mat','a_ini1','phi_ini1','gammaZ_ini1','a_ini2','phi_ini2','gammaZ_ini2','w_ini');
% vars = whos;
% for k = 1:length(vars)
%     name = vars(k).name;
%     val = eval(name);
%     disp([name ' = ' mat2str(val, 17) ';'])
% end

a_ini1 = [2.57894736842105e-06 1.63157894736842e-06 5.73684210526316e-06 4.1578947368421e-06 2.89473684210526e-06 6.05263157894737e-06 5.10526315789474e-06 3.21052631578947e-06 2.26315789473684e-06 6.68421052631579e-06 5.42105263157895e-06 1e-06 1.31578947368421e-06 1.94736842105263e-06 6.36842105263158e-06 4.78947368421053e-06 3.52631578947368e-06 7e-06 3.84210526315789e-06 4.47368421052632e-06];
a_ini2 = [6.47368421052632e-06 6.68421052631579e-06 8.36842105263158e-06 7.73684210526316e-06 8.57894736842105e-06 7.31578947368421e-06 5.21052631578947e-06 6.89473684210526e-06 8.1578947368421e-06 9e-06 5.42105263157895e-06 6.26315789473684e-06 8.78947368421053e-06 5.63157894736842e-06 5e-06 7.52631578947368e-06 7.10526315789474e-06 5.8421052631579e-06 7.94736842105263e-06 6.05263157894737e-06];
gammaZ_ini1 = [0.0418421052631579 0.0392105263157895 0.0155263157894737 0.055 0.0313157894736842 0.0497368421052632 0.00763157894736842 0.005 0.0207894736842105 0.0181578947368421 0.0444736842105263 0.0523684210526316 0.023421052631579 0.0365789473684211 0.0339473684210526 0.0260526315789474 0.0286842105263158 0.0102631578947368 0.0128947368421053 0.0471052631578947];
gammaZ_ini2 = [0.0523684210526316 0.0339473684210526 0.0418421052631579 0.005 0.0155263157894737 0.0365789473684211 0.0207894736842105 0.0102631578947368 0.055 0.0313157894736842 0.0128947368421053 0.0471052631578947 0.0181578947368421 0.023421052631579 0.0260526315789474 0.0286842105263158 0.00763157894736842 0.0444736842105263 0.0392105263157895 0.0497368421052632];
phi_ini1 = [0.173684210526316 0.3 0.0157894736842105 0.205263157894737 0.221052631578947 0.142105263157895 0.189473684210526 0.268421052631579 0.110526315789474 0.126315789473684 0 0.0631578947368421 0.252631578947368 0.157894736842105 0.0947368421052631 0.0789473684210526 0.236842105263158 0.0315789473684211 0.0473684210526316 0.284210526315789];
phi_ini2 = [0.734210526315789 0.623684210526316 0.607894736842105 0.686842105263158 0.528947368421053 0.75 0.45 0.576315789473684 0.718421052631579 0.639473684210526 0.497368421052632 0.592105263157895 0.702631578947368 0.481578947368421 0.560526315789474 0.544736842105263 0.655263157894737 0.465789473684211 0.513157894736842 0.671052631578947];
w_ini = [0.105263157894737 0.315789473684211 0 0.789473684210526 1 0.736842105263158 0.631578947368421 0.421052631578947 0.684210526315789 0.368421052631579 0.526315789473684 0.947368421052632 0.894736842105263 0.157894736842105 0.473684210526316 0.210526315789474 0.263157894736842 0.842105263157895 0.578947368421053 0.0526315789473684];


% a_ini=[4.42105263157895e-06,3.21052631578947e-06,1.77368421052632e-05,2.37894736842105e-05,2.50000000000000e-05,1.16842105263158e-05,1.41052631578947e-05,6.84210526315789e-06,8.05263157894737e-06,2.13684210526316e-05,2.01578947368421e-05,5.63157894736842e-06,1.89473684210526e-05,9.26315789473684e-06,1.28947368421053e-05,2.00000000000000e-06,1.04736842105263e-05,1.65263157894737e-05,1.53157894736842e-05,2.25789473684211e-05];
% phi_ini=[0.518421052631579,0.445789473684211,0.482105263157895,0.0100000000000000,0.409473684210526,0.627368421052632,0.554736842105263,0.591052631578947,0.118947368421053,0.336842105263158,0.264210526315790,0.191578947368421,0.700000000000000,0.155263157894737,0.0826315789473684,0.227894736842105,0.300526315789474,0.663684210526316,0.0463157894736842,0.373157894736842];
% gammaZ_ini=[0.150000000000000,0.0480526315789474,0.0323684210526316,0.126473684210526,0.0245263157894737,0.142157894736842,0.118631578947368,0.0872631578947369,0.0794210526315790,0.110789473684211,0.0558947368421053,0.0637368421052632,0.134315789473684,0.0951052631578947,0.0402105263157895,0.102947368421053,0.00100000000000000,0.0715789473684211,0.00884210526315789,0.0166842105263158];

options = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolX',1e-100);

for ii=1:length(a_ini1)
    % xinit=[1e-6 0.2 0.05];
    xinit=[
        a_ini1(ii), phi_ini1(ii), gammaZ_ini1(ii), ...
        a_ini2(ii), phi_ini2(ii), gammaZ_ini2(ii), ...
        w_ini(ii),
        ];
    [xmin,fvalmin,exitflag]=fminsearch(@(x) sfm2_F_myfun_SFMTroisParam(x,connu),xinit,options);
    % res=[xmin(1).*1e6 xmin(2) xmin(3)]
    if xmin(3)<0 || xmin(6)<0
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

%%%%%%%%%%%%%%


% xinit=[1e-6 0.6 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[1e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[1e-6 0.6 0.1];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% 
% xinit=[5e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% 
% xinit=[5e-6 0.6 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[5e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[5e-6 0.6 0.1];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% 
% xinit=[8e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% 
% xinit=[8e-6 0.6 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[8e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[8e-6 0.6 0.1];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% 
% xinit=[10e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% 
% xinit=[10e-6 0.6 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[10e-6 0.2 0.05];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xinit=[10e-6 0.6 0.1];
% [xmin,fvalmin,exitflag]=fminsearch(@(x) F_myfun_SFMTroisParam(x,connu),xinit,options);
% if xmin(3)<0
%    fvalmin=50;
% end
% if fvalmin<fvalminbest
%     xminbest=xmin;
%     fvalminbest=fvalmin;
% end
% 
% xminbest(4)=fvalminbest;

end