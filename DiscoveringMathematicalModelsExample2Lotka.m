%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 20;
dt  = 0.01;
t   = dt:dt:tf;

%% number of variables and coefficients
n = 2;
r = 14;

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 0;

%% state initialization
x        = [8;8];
xbar     = x;%zeros(n,1);
y        = x;zeros(n,1);
thetabar = zeros(r,1);
 
%% true parameters
alpha = 1;
beta = 0.2;
delta = 0.1;
gamma = 0.2;

%% initial control inputs
%u     = [10 1]';

%% for plotting
%uArray          = [];
xArray          = [];
xbarArray       = [];
yArray          = [];
thetabarArray   = [];

%% Initialization for estimator

lambdav = 0.995;
lambdat = 0.999;
Rx      = 1*eye(n);
Rt      = 1*eye(n);
Px      = 0.1*eye(n);
Pt      = 0.1*eye(r);
Gamma   = 1*zeros(n,r);

%% simulation
for i=1:(tf/dt)

    %uArray         = [uArray u];
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];    
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 
    
    x = A*x+dt*[alpha*x(1)-beta*x(1)*x(2);delta*x(1)*x(2)-gamma*x(2)];
    y = C*x+dt*R^2*randn(n,1);

    Phi = [y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3 zeros(7,1)';
           zeros(7,1)' y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3];
    
    % Estimation using adaptive observer
    Kx = Px*C'*inv(C*Px*C'+Rx);
    Kt = Pt*Gamma'*C'*inv(C*Gamma*Pt*Gamma'*C'+Rt);
    Gamma = (eye(n)-Kx*C)*Gamma;

    xbar = xbar+(Kx+Gamma*Kt)*(y-C*xbar);
    thetabar = thetabar-Kt*(y-C*xbar);

    xbar = A*xbar+Phi*thetabar;

    thetabar = thetabar;
    Px = (1/lambdav)*eye(n)*(eye(n)-Kx*C)*Px*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*C*Gamma)*Pt;
    Gamma = eye(n)*Gamma-Phi;

end

figure(1)
plot(t,yArray(1,:),'-','LineWidth',16);
hold on;
plot(t,xbarArray(1,:),':','LineWidth',16)
hold on;
plot(t,yArray(2,:),'-','LineWidth',16);
hold on;
plot(t,xbarArray(2,:),':','LineWidth',16)
set(gca,'color','white','LineWidth',3,'FontSize',56)
grid on;
grid minor;
xlabel('t (s)')

figure(2)
subplot(4,2,1)
plot(t,alpha*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true','estimated')
grid on;
grid minor;
ylim([alpha-0.2 alpha+0.2])
ylabel('\alpha','FontSize',72)
subplot(4,2,2)
plot(t,alpha*ones(1,length(t))-thetabarArray(1,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \alpha','FontSize',72)
subplot(4,2,3)
plot(t,beta*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(5,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0 0.8])
ylabel('\beta','FontSize',72)
subplot(4,2,4)
plot(t,beta*ones(1,length(t))+thetabarArray(5,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \beta','FontSize',72)
subplot(4,2,5)
plot(t,delta*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(12,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0 0.2])
ylabel('\delta','FontSize',72)
subplot(4,2,6)
plot(t,delta*ones(1,length(t))-thetabarArray(12,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \delta','FontSize',72)
subplot(4,2,7)
plot(t,gamma*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(9,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0 0.8])
ylabel('\gamma','FontSize',72)
xlabel('t (s)')
subplot(4,2,8)
plot(t,gamma*ones(1,length(t))+thetabarArray(9,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \gamma','FontSize',72)
xlabel('t (s)')

Coeff = round([(1/dt)*thetabar(1:(r/n),end)'; (1/dt)*thetabar((r/n)+1:2*(r/n),end)'; (1/dt)*thetabar(2*(r/n)+1:r,end)'],2)
