clear
clc
close all

%% Including Figure 6D 6E 6F and SI Figure 9 A-H

gamma_BA_test=0:0.002:0.2;
gamma_BD_test=0:0.002:0.2;

r_GA=1;
r_GD=1;
r_BA=1;
r_BD=1;

K_GA=0.2;
K_GD=0.2;
K_BA_test=0.2:0.02:2;
K_BD_test=0.2:0.02:2;

r_GP=1;
r_BC=1;

K_GP=0.2;
K_BC=0.2;

eta_test=0:0.01:1;
tspan=0:0.05:96;

cmap(1:128,1)=0.85:0.15/127:1;
cmap(1:128,2)=0.33:0.67/127:1;
cmap(1:128,3)=0.10:0.90/127:1;
cmap(129:256,1)=1:-1.00/127:0.00;
cmap(129:256,2)=1:-0.55/127:0.45;
cmap(129:256,3)=1:-0.25/127:0.75;
cmap=cmap(end:-1:1,:);

P_test=0.0:0.02:2;
C_test=2:-0.02:0;



%% eta_test Figure 6D

for i=1:size(P_test,2)
    gamma_GA=0.1;
    gamma_GD=0.1;
    gamma_BA=0.001;
    gamma_BD=0.15;
    K_BA=2;    
    K_BD=0.5;
    eta=0.25;

    InitP=P_test(i);
    InitC=C_test(i);

    for j=1:size(eta_test,2)

        eta=eta_test(j);

        c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

        InitCo=[0.05,0.05,InitP,InitC,0,0]';
        InitG=[0.05,0,InitP,InitC,0,0]';
        InitB=[0,0.05,InitP,InitC,0,0]';

        [~,TestCo] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCo);
        [~,TestG] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitG);
        [~,TestB] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitB);

        InteractionBtoG(j,i)=sum(TestCo(:,1)-TestG(:,1))./1000;
        InteractionGtoB(j,i)=sum(TestCo(:,2)-TestB(:,2))./1000;
        Consen_Gal(j,i)=TestCo(361,5);
        Consen_diGal(j,i)=TestCo(361,6);

    end

    DisplayBar(i, size(P_test,2));
end


[X,Y] = meshgrid(P_test./2,eta_test);

figure('Position',[100,100,1600,400],'Color',[1,1,1])
subplot(1,3,1)
P1=pcolor(X,Y,InteractionBtoG);
P1.LineStyle='none';
P1.FaceColor = 'interp';
axis([min(P_test./2),max(P_test./2),min(eta_test),max(eta_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('PGA Percentage')
ylabel('\eta')
title('Interaction Bi to Ga')


subplot(1,3,2)
P2=pcolor(X,Y,InteractionGtoB);
P2.LineStyle='none';
P2.FaceColor = 'interp';
axis([min(P_test./2),max(P_test./2),min(eta_test),max(eta_test)])

box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('PGA Percentage')
ylabel('\eta')
title('Interaction Ga to Bi')

subplot(1,3,3) % figure 6D
P3=pcolor(X,Y,log2(Consen_Gal./Consen_diGal));
P3.LineStyle='none';
P3.FaceColor = 'interp';
axis([min(P_test./2),max(P_test./2),min(eta_test),max(eta_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-4 4])
colormap(cmap)
colorbar
xlabel('PGA Percentage')
ylabel('\eta')
title('log_2(Gal / diGal) (18h)')
set(gcf,'PaperType','A2')

% exportgraphics(gcf,'eta.pdf','ContentType','Vector')

Consen_Gal=[];
Consen_diGal=[];
for i=1:size(P_test,2)
    gamma_GA=0.1;
    gamma_GD=0.1;
    gamma_BA=0.001;
    gamma_BD=0.15;
    K_BA=2;
    K_BD=0.5;
    eta=0.25;

    InitP=P_test(i);
    InitC=C_test(i);

    c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

    InitCo=[0.05,0.05,InitP,InitC,0,0]';
    InitG=[0.05,0,InitP,InitC,0,0]';
    InitB=[0,0.05,InitP,InitC,0,0]';

    [~,TestCo] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCo);
    [~,TestG] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitG);
    [~,TestB] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitB);


    Consen_Gal(i)=TestCo(361,5);
    Consen_diGal(i)=TestCo(361,6);
end

figure('Position',[100,100,800,500],'Color',[1,1,1])
hold on
plot(0:0.01:1,Consen_Gal,'LineWidth',1.5);
plot(0:0.01:1,Consen_diGal,'LineWidth',1.5);

axis([min(P_test./2),max(P_test./2),0,1.2])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
xlabel('PGA Percentage')
ylabel('Consentration (a.u.)')
legend('GalA','diGalA','Location','northwest')
title('Interaction Bi to Ga')

% exportgraphics(gcf,'GalAdiGalA.pdf','ContentType','Vector')

%% gamma_BA vs K_BA PGA SI Figure 9 A and E

InteractionBtoG=[];
InteractionGtoB=[];
Consen_Gal=[];
Consen_diGal=[];

for i=1:size(gamma_BA_test,2)
    gamma_GA=0.1;
    gamma_GD=0.1;
    gamma_BD=0.15;
    K_BD=0.5;
    eta=0.25;

    InitP=P_test(end);
    InitC=C_test(end);

    gamma_BA=gamma_BA_test(i);

    for j=1:size(K_BA_test,2)

        K_BA=K_BA_test(j);

        c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

        InitCo=[0.05,0.05,InitP,InitC,0,0]';
        InitG=[0.05,0,InitP,InitC,0,0]';
        InitB=[0,0.05,InitP,InitC,0,0]';

        [~,TestCo] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCo);
        [~,TestG] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitG);
        [~,TestB] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitB);

        InteractionBtoG(j,i)=sum(TestCo(:,1)-TestG(:,1))./1000;
        InteractionGtoB(j,i)=sum(TestCo(:,2)-TestB(:,2))./1000;
        Consen_Gal(j,i)=TestCo(361,5);
        Consen_diGal(j,i)=TestCo(361,6);

    end

    DisplayBar(i,size(gamma_BA_test,2));
end
InteractionBtoG_PGA_BA=InteractionBtoG;
InteractionGtoB_PGA_BA=InteractionGtoB;

[X,Y] = meshgrid(gamma_BA_test,K_BA_test);

figure('Position',[100,100,1600,400],'Color',[1,1,1])
subplot(1,3,1)
hold on
P1=contourf(X,Y,InteractionBtoG);
axis([min(gamma_BA_test),max(gamma_BA_test),min(K_BA_test),max(K_BA_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('Interaction Bi to Ga')


subplot(1,3,2)
P2=contourf(X,Y,InteractionGtoB);
axis([min(gamma_BA_test),max(gamma_BA_test),min(K_BA_test),max(K_BA_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('Interaction Ga to Bi')

subplot(1,3,3)
P3=pcolor(X,Y,log2(Consen_Gal./Consen_diGal));
P3.LineStyle='none';
P3.FaceColor = 'interp';
axis([min(gamma_BA_test),max(gamma_BA_test),min(K_BA_test),max(K_BA_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-1 1])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('log_2(Gal / diGal) (18h)')

set(gcf,'PaperType','A2')

% exportgraphics(gcf,'PGA_BA_test.pdf','ContentType','Vector')

%% gamma_BD vs K_BD PGA  SI Figure 9 B and F

InteractionBtoG=[];
InteractionGtoB=[];
Consen_Gal=[];
Consen_diGal=[];

for i=1:size(gamma_BD_test,2)
    gamma_GA=0.1;
    gamma_GD=0.1;
    gamma_BA=0.001;
    gamma_BD=0.15;
    K_BA=2;    
    K_BD=0.5;
    eta=0.25;

    InitP=P_test(end);
    InitC=C_test(end);

    gamma_BD=gamma_BD_test(i);

    for j=1:size(K_BD_test,2)
        
        K_BD=K_BD_test(j);

        c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

        InitCo=[0.05,0.05,InitP,InitC,0,0]';
        InitG=[0.05,0,InitP,InitC,0,0]';
        InitB=[0,0.05,InitP,InitC,0,0]';

        [~,TestCo] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCo);
        [~,TestG] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitG);
        [~,TestB] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitB);

        InteractionBtoG(j,i)=sum(TestCo(:,1)-TestG(:,1))./1000;
        InteractionGtoB(j,i)=sum(TestCo(:,2)-TestB(:,2))./1000;
        Consen_Gal(j,i)=TestCo(361,5);
        Consen_diGal(j,i)=TestCo(361,6);

    end

    DisplayBar(i,size(gamma_BD_test,2));
end
InteractionBtoG_PGA_BD=InteractionBtoG;
InteractionGtoB_PGA_BD=InteractionGtoB;

[X,Y] = meshgrid(gamma_BD_test,K_BD_test);

figure('Position',[100,100,1600,400],'Color',[1,1,1])
subplot(1,3,1)
P1=contourf(X,Y,InteractionBtoG);

axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('gamma_{BD}')
ylabel('K_{BD}')
title('Interaction Bi to Ga')


subplot(1,3,2)
P2=contourf(X,Y,InteractionGtoB);
% P2.LineStyle='none';
% P2.FaceColor = 'interp';
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('gamma_{BD}')
ylabel('K_{BD}')
title('Interaction Ga to Bi')

subplot(1,3,3)
P3=pcolor(X,Y,log2(Consen_Gal./Consen_diGal));
P3.LineStyle='none';
P3.FaceColor = 'interp';
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-1 1])
colormap(cmap)
colorbar
xlabel('gamma_{BD}')
ylabel('K_{BD}')
title('log_2(Gal / diGal) (18h)')

set(gcf,'PaperType','A2')

% exportgraphics(gcf,'PGA_BD_test.pdf','ContentType','Vector')

%% gamma_BA vs K_BA PC  SI Figure 9 C and G

InteractionBtoG=[];
InteractionGtoB=[];
Consen_Gal=[];
Consen_diGal=[];


for i=1:size(gamma_BA_test,2)
    gamma_GA=0.1;
    gamma_GD=0.1;
    gamma_BA=0.001;
    gamma_BD=0.15;
    K_BA=2;    
    K_BD=0.5;
    eta=0.25;

    InitP=P_test(31);
    InitC=C_test(31);

    gamma_BA=gamma_BA_test(i);

    for j=1:size(K_BA_test,2)

        K_BA=K_BA_test(j);

        c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

        InitCo=[0.05,0.05,InitP,InitC,0,0]';
        InitG=[0.05,0,InitP,InitC,0,0]';
        InitB=[0,0.05,InitP,InitC,0,0]';

        [~,TestCo] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCo);
        [~,TestG] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitG);
        [~,TestB] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitB);

        InteractionBtoG(j,i)=sum(TestCo(:,1)-TestG(:,1))./1000;
        InteractionGtoB(j,i)=sum(TestCo(:,2)-TestB(:,2))./1000;
        Consen_Gal(j,i)=TestCo(361,5);
        Consen_diGal(j,i)=TestCo(361,6);

    end

    DisplayBar(i,size(gamma_BA_test,2));
end
InteractionBtoG_PC_BA=InteractionBtoG;
InteractionGtoB_PC_BA=InteractionGtoB;
[X,Y] = meshgrid(gamma_BA_test,K_BA_test);

figure('Position',[100,100,1600,400],'Color',[1,1,1])
subplot(1,3,1)
hold on
P1=contourf(X,Y,InteractionBtoG);
axis([min(gamma_BA_test),max(gamma_BA_test),min(K_BA_test),max(K_BA_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('Interaction Bi to Ga')


subplot(1,3,2)
P2=contourf(X,Y,InteractionGtoB);
axis([min(gamma_BA_test),max(gamma_BA_test),min(K_BA_test),max(K_BA_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('Interaction Ga to Bi')

subplot(1,3,3)
P3=pcolor(X,Y,log2(Consen_Gal./Consen_diGal));
P3.LineStyle='none';
P3.FaceColor = 'interp';
axis([min(gamma_BA_test),max(gamma_BA_test),min(K_BA_test),max(K_BA_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-1 1])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('log_2(Gal / diGal) (18h)')

set(gcf,'PaperType','A2')
% exportgraphics(gcf,'PC_BA_test.pdf','ContentType','Vector')

%% gamma_BD vs K_BD PC  SI Figure 9 D and H

InteractionBtoG=[];
InteractionGtoB=[];
Consen_Gal=[];
Consen_diGal=[];

for i=1:size(gamma_BD_test,2)
    gamma_GA=0.1;
    gamma_GD=0.1;
    gamma_BA=0.001;
    gamma_BD=0.15;
    K_BA=2;    
    K_BD=0.5;
    eta=0.25;

    InitP=P_test(31);
    InitC=C_test(31);

    gamma_BD=gamma_BD_test(i);

    for j=1:size(K_BD_test,2)
        
        K_BD=K_BD_test(j);

        c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

        InitCo=[0.05,0.05,InitP,InitC,0,0]';
        InitG=[0.05,0,InitP,InitC,0,0]';
        InitB=[0,0.05,InitP,InitC,0,0]';

        [~,TestCo] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCo);
        [~,TestG] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitG);
        [~,TestB] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitB);

        InteractionBtoG(j,i)=sum(TestCo(:,1)-TestG(:,1))./1000;
        InteractionGtoB(j,i)=sum(TestCo(:,2)-TestB(:,2))./1000;
        Consen_Gal(j,i)=TestCo(361,5);
        Consen_diGal(j,i)=TestCo(361,6);

    end

    DisplayBar(i,size(gamma_BD_test,2));
end
InteractionBtoG_PC_BD=InteractionBtoG;
InteractionGtoB_PC_BD=InteractionGtoB;
[X,Y] = meshgrid(gamma_BD_test,K_BD_test);

figure('Position',[100,100,1600,400],'Color',[1,1,1])
subplot(1,3,1)
P1=contourf(X,Y,InteractionBtoG);
% P1.LineStyle='none';
% P1.FaceColor = 'interp';
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('gamma_{BD}')
ylabel('K_{BD}')
title('Interaction Bi to Ga')


subplot(1,3,2)
P2=contourf(X,Y,InteractionGtoB);
% P2.LineStyle='none';
% P2.FaceColor = 'interp';
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-2 2])
colormap(cmap)
colorbar
xlabel('gamma_{BD}')
ylabel('K_{BD}')
title('Interaction Ga to Bi')

subplot(1,3,3)
P3=pcolor(X,Y,log2(Consen_Gal./Consen_diGal));
P3.LineStyle='none';
P3.FaceColor = 'interp';
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-1 1])
colormap(cmap)
colorbar
xlabel('gamma_{BD}')
ylabel('K_{BD}')
title('log_2(Gal / diGal) (18h)')

set(gcf,'PaperType','A2')
% exportgraphics(gcf,'PC_BD_test.pdf','ContentType','Vector')





%% Figure 6E

figure('Position',[100,100,600,500],'Color',[1,1,1])
contourf(X,Y,InteractionGtoB_PC_BA./InteractionGtoB_PGA_BA)
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-4 4])
colormap(cmap)
colorbar
xlabel('\gamma_{BA}')
ylabel('K_{BA}')
title('PC v.s. PGA')

set(gcf,'PaperType','A2')
% exportgraphics(gcf,'GtoB_PCvsPGA_BA_test.pdf','ContentType','Vector')



%% Figure 6F

figure('Position',[100,100,600,500],'Color',[1,1,1])
contourf(X,Y,InteractionGtoB_PC_BD./InteractionGtoB_PGA_BD)
axis([min(gamma_BD_test),max(gamma_BD_test),min(K_BD_test),max(K_BD_test)])
box on
set(gca,'Fontsize',14,'LineWidth',2)
set(gca,'Layer','top')
set(gca,'clim',[-4 4])
colormap(cmap)
colorbar
xlabel('\gamma_{BD}')
ylabel('K_{BD}')
title('PC v.s. PGA')

set(gcf,'PaperType','A2')
% exportgraphics(gcf,'GtoB_PCvsPGA_BD_test.pdf','ContentType','Vector')

