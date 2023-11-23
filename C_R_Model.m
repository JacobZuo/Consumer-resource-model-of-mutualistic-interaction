function dydt = C_R_Model(t,y,c)
x_G=y(1);
x_B=y(2);

P=y(3);
C=y(4);

A=y(5);
D=y(6);

gamma_GA=c(1);
gamma_GD=c(2);
gamma_BA=c(3);
gamma_BD=c(4);

r_GA=c(5);
r_GD=c(6);
r_BA=c(7);
r_BD=c(8);

K_GA=c(9);
K_GD=c(10);
K_BA=c(11);
K_BD=c(12);

r_GP=c(13);
r_BC=c(14);
K_GP=c(15);
K_BC=c(16);

r_P=0.005;
r_C=0.005;
K_P=20;
K_C=20;
r_D=0.005;
K_D=20;


eta=c(17);

dx_G=(gamma_GA.*r_GA.*A./(K_GA+A)+gamma_GD.*r_GD.*D./(K_GD+D)).*x_G;
dx_B=(gamma_BA.*r_BA.*A./(K_BA+A)+gamma_BD.*r_BD.*D./(K_BD+D)).*x_B;

dP=-r_GP.*P./(K_GP+P).*x_G-r_P.*P./(K_P+P)+r_BC.*C./(K_BC+C).*x_B+r_C.*C./(K_C+C);
dC=-r_BC.*C./(K_BC+C).*x_B-r_C.*C./(K_C+C);

dA=-gamma_GA.*r_GA.*A./(K_GA+A).*x_G-gamma_BA.*r_BA.*A./(K_BA+A).*x_B+(1./((eta./(P./2)).^2+1)).*(r_GP.*P./(K_GP+P).*x_G+r_P.*P./(K_P+P))+r_D.*D./(K_D+D);
dD=-gamma_GD.*r_GD.*D./(K_GD+D).*x_G-gamma_BD.*r_BD.*D./(K_BD+D).*x_B+(1-(1./((eta./(P./2)).^2+1))).*(r_GP.*P./(K_GP+P).*x_G+r_P.*P./(K_P+P))-r_D.*D./(K_D+D);



dydt=[dx_G,dx_B,dP,dC,dA,dD]';


end