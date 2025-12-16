function[Params]=sfm2_inversion_BSC_SFM_Neldermead_sansLog_Fc(freq,data,nb_xini)

% data=load('bsc_freq.mat');

BSC=data;

soundspeed=1540;
kk=2*pi*freq/soundspeed;
connu=[BSC kk];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% '_______________________________________'
% 'Neldermead - SFM Three parameters - SANS LOG : [a phi gammaZ CdB]'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


[x,fvalmin,exitflag]=sfm2_F_InversionSFMTroisParam_NelderMead_Fc(connu,nb_xini);
Vs=(4/3)*pi*(x(1))^3;
% CdB=10*log10((x(3)^2*x(2)/Vs)*(1e-3)^3); % is complex in a few cases
CdB=0;

% BSC_opt=sfm2_F_FluidSFMMonodisperse(kk,x(1),x(2),x(3)); TODO

% R_Neldermead = 1-(sum((BSC-BSC_opt).^2)/sum((BSC-mean(BSC)).^2)); TODO
R_Neldermead=0;
% % 
% figure
% plot(freq,BSC,'r')
% hold on
% plot(freq,BSC_opt,'b')

a1=x(1);
phi1=x(2);
gammaZ1=x(3);
a2=x(4);
phi2=x(5);
gammaZ2=x(6);
w=x(7);

clear x
clear BSC_opt

Params=[a1 phi1 gammaZ1 a2 phi2 gammaZ2 w CdB R_Neldermead fvalmin exitflag];


end
