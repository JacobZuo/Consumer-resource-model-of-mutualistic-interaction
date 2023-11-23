clear
clc
close all

%% Including Figure 6G

gamma_GA=0.1;
gamma_GD=0.1;
gamma_BA=0.001;
gamma_BD=0.05;

K_GA=0.2;
K_GD=0.2;
K_BA=2;
K_BD=2;

K_GP=0.2;
K_BC=0.2;

eta=0.25;
tspan=0:0.05:96;


% The following parameters set to 1 would not affect the simulation
r_GA=1;
r_GD=1;
r_BA=1;
r_BD=1;
r_GP=1;
r_BC=1;


%%
c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

InitCoOnP=[0.05,0.05,2,0,0,0]';
InitGOnP=[0.05,0,2,0,0,0]';
InitBOnP=[0,0.05,2,0,0,0]';

InitCoOnC=[0.05,0.05,0.5,1.5,0,0]';
InitGOnC=[0.05,0,0.6,1.4,0,0]';
InitBOnC=[0,0.05,0.6,1.4,0,0]';

[~,TestCoOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCoOnP);
[~,TestGOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitGOnP);
[~,TestBOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitBOnP);

[~,TestCoOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCoOnC);
[~,TestGOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitGOnC);
[~,TestBOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitBOnC);


AUCG_OnP=sum(TestGOnP(:,1))./1000;
AUCG_CoOnP=sum(TestCoOnP(:,1))./1000;
AUCB_OnP=sum(TestBOnP(:,2))./1000;
AUCB_CoOnP=sum(TestCoOnP(:,2))./1000;

AUCG_OnC=sum(TestGOnC(:,1))./1000;
AUCG_CoOnC=sum(TestCoOnC(:,1))./1000;
AUCB_OnC=sum(TestBOnC(:,2))./1000;
AUCB_CoOnC=sum(TestCoOnC(:,2))./1000;


figure('Position',[100,100,600,400],'Color',[1,1,1])
hold on
plot(tspan,TestBOnC(:,2),'Color',[0.6 0.75 1],'LineWidth',2)
plot(tspan,TestCoOnC(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)

plot(tspan,TestBOnP(:,2),':','Color',[0.4 0.55 0.8],'LineWidth',2)
plot(tspan,TestCoOnP(:,2),':','Color',[0 0.3470 0.541],'LineWidth',2)
% set(gca,'YScale','log')
axis([0 96 0.0 0.2])
xlabel('time (h)')
ylabel('Biomass (a.u. log)')
box on
legend('Bi in PC','Bi(BG) in PC','Bi in PGA','Bi(BG) in PGA','Location','Northwest')
title('PGA')
set(gca,'Fontsize',14,'LineWidth',2)
set(gcf,'PaperType','A2')
% print(['Figure 6G_1.pdf'],'-dpdf','-r300')

%%
gamma_BD=0.15;
K_BD=0.5;

c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

InitCoOnP=[0.05,0.05,2,0,0,0]';
InitGOnP=[0.05,0,2,0,0,0]';
InitBOnP=[0,0.05,2,0,0,0]';

InitCoOnC=[0.05,0.05,0.5,1.5,0,0]';
InitGOnC=[0.05,0,0.6,1.4,0,0]';
InitBOnC=[0,0.05,0.6,1.4,0,0]';

[~,TestCoOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCoOnP);
[~,TestGOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitGOnP);
[~,TestBOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitBOnP);

[~,TestCoOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCoOnC);
[~,TestGOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitGOnC);
[~,TestBOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitBOnC);


AUCG_OnP=sum(TestGOnP(:,1))./1000;
AUCG_CoOnP=sum(TestCoOnP(:,1))./1000;
AUCB_OnP=sum(TestBOnP(:,2))./1000;
AUCB_CoOnP=sum(TestCoOnP(:,2))./1000;

AUCG_OnC=sum(TestGOnC(:,1))./1000;
AUCG_CoOnC=sum(TestCoOnC(:,1))./1000;
AUCB_OnC=sum(TestBOnC(:,2))./1000;
AUCB_CoOnC=sum(TestCoOnC(:,2))./1000;


figure('Position',[100,100,600,400],'Color',[1,1,1])
hold on
plot(tspan,TestBOnC(:,2),'Color',[0.6 0.75 1],'LineWidth',2)
plot(tspan,TestCoOnC(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)

plot(tspan,TestBOnP(:,2),':','Color',[0.4 0.55 0.8],'LineWidth',2)
plot(tspan,TestCoOnP(:,2),':','Color',[0 0.3470 0.541],'LineWidth',2)
% set(gca,'YScale','log')
axis([0 96 0.0 0.5])
xlabel('time (h)')
ylabel('Biomass (a.u. log)')
box on
legend('Bi in PC','Bi(BG) in PC','Bi in PGA','Bi(BG) in PGA','Location','Northwest')
title('PGA')
set(gca,'Fontsize',14,'LineWidth',2)
set(gcf,'PaperType','A2')
% print(['Figure 6G_2.pdf'],'-dpdf','-r300')




%%
gamma_BD=0.2;
K_BD=0.2;
c=[gamma_GA,gamma_GD,gamma_BA,gamma_BD,r_GA,r_GD,r_BA,r_BD,K_GA,K_GD,K_BA,K_BD,r_GP,r_BC,K_GP,K_BC,eta];

InitCoOnP=[0.05,0.05,2,0,0,0]';
InitGOnP=[0.05,0,2,0,0,0]';
InitBOnP=[0,0.05,2,0,0,0]';

InitCoOnC=[0.05,0.05,0.5,1.5,0,0]';
InitGOnC=[0.05,0,0.6,1.4,0,0]';
InitBOnC=[0,0.05,0.6,1.4,0,0]';

[~,TestCoOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCoOnP);
[~,TestGOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitGOnP);
[~,TestBOnP] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitBOnP);

[~,TestCoOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitCoOnC);
[~,TestGOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitGOnC);
[~,TestBOnC] = ode45(@(t,y) C_R_Model(t,y,c),tspan,InitBOnC);


AUCG_OnP=sum(TestGOnP(:,1))./1000;
AUCG_CoOnP=sum(TestCoOnP(:,1))./1000;
AUCB_OnP=sum(TestBOnP(:,2))./1000;
AUCB_CoOnP=sum(TestCoOnP(:,2))./1000;

AUCG_OnC=sum(TestGOnC(:,1))./1000;
AUCG_CoOnC=sum(TestCoOnC(:,1))./1000;
AUCB_OnC=sum(TestBOnC(:,2))./1000;
AUCB_CoOnC=sum(TestCoOnC(:,2))./1000;


figure('Position',[100,100,600,400],'Color',[1,1,1])
hold on
plot(tspan,TestBOnC(:,2),'Color',[0.6 0.75 1],'LineWidth',2)
plot(tspan,TestCoOnC(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)

plot(tspan,TestBOnP(:,2),':','Color',[0.4 0.55 0.8],'LineWidth',2)
plot(tspan,TestCoOnP(:,2),':','Color',[0 0.3470 0.541],'LineWidth',2)
% set(gca,'YScale','log')
axis([0 96 0.0 1.2])
xlabel('time (h)')
ylabel('Biomass (a.u. log)')
box on
legend('Bi in PC','Bi(BG) in PC','Bi in PGA','Bi(BG) in PGA','Location','Northwest')
title('PGA')
set(gca,'Fontsize',14,'LineWidth',2)
set(gcf,'PaperType','A2')
% print(['Figure 6G_3.pdf'],'-dpdf','-r300')


%%

% figure('Position',[100,100,1200,800],'Color',[1,1,1])
% 
% subplot(2,3,1)
% hold on
% plot(tspan,TestGOnC(:,1),'Color',[0.85 0.325 0.098],'LineWidth',2)
% plot(tspan,TestBOnC(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)
% plot(tspan,TestCoOnC(:,1)+TestCoOnC(:,2),'Color',[0.466 0.674 0.188],'LineWidth',2)
% axis([0 96 0 3])
% xlabel('time (h)')
% ylabel('Biomass (a.u.)')
% box on
% legend('Ga','Bi','BG','Location','Northwest')
% set(gca,'Fontsize',14,'LineWidth',2)
% 
% title('PC')
% 
% subplot(2,3,2)
% hold on
% plot(tspan,TestGOnC(:,1),'Color',[1 0.75 0.5],'LineWidth',2)
% plot(tspan,TestCoOnC(:,1),'Color',[0.85 0.325 0.098],'LineWidth',2)
% % set(gca,'YScale','log')
% axis([0 96 0.0 3])
% xlabel('time (h)')
% ylabel('Biomass (a.u. log)')
% box on
% legend('Ga','Ga(BG)','Location','Northwest')
% title('PC')
% set(gca,'Fontsize',14,'LineWidth',2)
% 
% 
% subplot(2,3,3)
% hold on
% plot(tspan,TestBOnC(:,2),'Color',[0.6 0.75 1],'LineWidth',2)
% plot(tspan,TestCoOnC(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)
% % set(gca,'YScale','log')
% axis([0 96 0.0 0.5])
% xlabel('time (h)')
% ylabel('Biomass (a.u. log)')
% box on
% legend('Bi','Bi(BG)','Location','Northwest')
% title('PC')
% set(gca,'Fontsize',14,'LineWidth',2)
% 
% 
% subplot(2,3,4)
% 
% hold on
% plot(tspan,TestGOnP(:,1),'Color',[0.85 0.325 0.098],'LineWidth',2)
% plot(tspan,TestBOnP(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)
% plot(tspan,TestCoOnP(:,1)+TestCoOnP(:,2),'Color',[0.466 0.674 0.188],'LineWidth',2)
% axis([0 96 0 3])
% xlabel('time (h)')
% ylabel('Biomass (a.u.)')
% box on
% legend('Ga','Bi','BG','Location','Northwest')
% title('PGA')
% set(gca,'Fontsize',14,'LineWidth',2)
% 
% subplot(2,3,5)
% 
% hold on
% plot(tspan,TestGOnP(:,1),'Color',[1 0.75 0.5],'LineWidth',2)
% plot(tspan,TestCoOnP(:,1),'Color',[0.85 0.325 0.098],'LineWidth',2)
% % set(gca,'YScale','log')
% axis([0 96 0.0 3])
% xlabel('time (h)')
% ylabel('Biomass (a.u. log)')
% box on
% legend('Ga','Ga(BG)','Location','Northwest')
% title('PGA')
% set(gca,'Fontsize',14,'LineWidth',2)
% 
% 
% subplot(2,3,6)
% 
% hold on
% plot(tspan,TestBOnP(:,2),'Color',[0.6 0.75 1],'LineWidth',2)
% plot(tspan,TestCoOnP(:,2),'Color',[0 0.4470 0.741],'LineWidth',2)
% % set(gca,'YScale','log')
% axis([0 96 0.0 0.5])
% xlabel('time (h)')
% ylabel('Biomass (a.u. log)')
% box on
% legend('Bi','Bi(BG)','Location','Northwest')
% title('PGA')
% set(gca,'Fontsize',14,'LineWidth',2)
% set(gcf,'PaperType','A2')



