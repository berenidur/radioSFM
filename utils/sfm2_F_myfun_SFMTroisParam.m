function [F] = sfm2_F_myfun_SFMTroisParam(x,connu) % Structure Factor Model
% x(1) scatterer radius         1
% x(2) scatterer concentration  1
% x(3) impedance contrast       1
% x(4) scatterer radius         2
% x(5) scatterer concentration  2
% x(6) impedance contrast       2
% x(7) n to c scattering ratio  w

kk=connu(:,2);
w=x(7);

% Facteur de structure analytique (hard sphere 3D)
% Voir Annexe de l'article Franceschini & Guillermin JASA 2012
% References : Ashcroft Phys Rev 1966 et Wherteim PRL 1963
q1=2*kk*x(1);
integrale1=(3*x(2)*(x(2)/2 + 1)^2*(cos(2*q1) + 2*q1.*sin(2*q1) - 2*q1.^2.*cos(2*q1) - 1))./(4*q1.^4*(x(2) - 1)^4) - ((2*x(2) + 1)^2*(sin(2*q1) - 2*q1.*cos(2*q1)))./(8*q1.^3*(x(2) - 1)^4) - (x(2)*(2*x(2) + 1)^2*(sin(2*q1).*(q1.^(-2) - (3/2)*q1.^(-4)) - cos(2*q1).*((1/2)*q1.^(-1) - (3/2)*q1.^(-3) + (3/4)*q1.^(-5)) + (3/4)*q1.^(-5)))./(4*q1.*(x(2) - 1)^4);
S1=ones(size(q1,2),1)'./(1-4*pi*(x(2)/(4*pi*x(1)^3/3))*(2*x(1))^3*integrale1);                  % Facteur de structure  

q2=2*kk*x(4);
integrale2=(3*x(5)*(x(5)/2 + 1)^2*(cos(2*q2) + 2*q2.*sin(2*q2) - 2*q2.^2.*cos(2*q2) - 1))./(4*q2.^4*(x(5) - 1)^4) - ((2*x(5) + 1)^2*(sin(2*q2) - 2*q2.*cos(2*q2)))./(8*q2.^3*(x(5) - 1)^4) - (x(5)*(2*x(5) + 1)^2*(sin(2*q2).*(q2.^(-2) - (3/2)*q2.^(-4)) - cos(2*q2).*((1/2)*q2.^(-1) - (3/2)*q2.^(-3) + (3/4)*q2.^(-5)) + (3/4)*q2.^(-5)))./(4*q2.*(x(5) - 1)^4);
S2=ones(size(q2,2),1)'./(1-4*pi*(x(5)/(4*pi*x(4)^3/3))*(2*x(4))^3*integrale2);                  % Facteur de structure

% Modele de Facteur de Structure
% BSC(k)=m sigmab(k) S(k)
% BSC(k)=(phi/Vs) (k^4 Vs^2 gammaz^2/(4 pi^2)) FormFactor  S(k)
% BSC(k)= (phi k^4 Vs gammaz^2/(4 pi^2)) FormFactor  S(k)
FormFactor1=((3*(sin(2.*kk.*x(1))-2.*kk.*x(1).*cos(2.*kk.*x(1)))./(2.*kk.*x(1)).^3).^2);
FormFactor2=((3*(sin(2.*kk.*x(4))-2.*kk.*x(4).*cos(2.*kk.*x(4)))./(2.*kk.*x(4)).^3).^2);
% 

bsc1=(1/(4*pi^2))*(4/3)*pi*x(1)^3*x(2)*x(3)^2*kk.^4.*FormFactor1.*S1;
bsc2=(1/(4*pi^2))*(4/3)*pi*x(4)^3*x(5)*x(6)^2*kk.^4.*FormFactor2.*S2;
bscnc=w*bsc1+(1-w)*bsc2;

F=sum((connu(:,1)-bscnc).^2)/sum(connu(:,1).^2);


end

