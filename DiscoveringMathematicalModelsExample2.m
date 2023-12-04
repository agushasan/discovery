%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 4;
dt  = 0.001;
t   = dt:dt:tf;

%% number of variables and coefficients
n = 3;
r = 48;

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 1;

%% state initialization
x        = [-8;7;27];
xbar     = zeros(n,1);
y        = zeros(n,1);
thetabar = zeros(r,1);
 
%% true parameters
sigma   = 10;
rho     = 28;
beta    = 3;

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
    
    x = A*x+dt*[10*(x(2)-x(1));x(1)*(28-x(3))-x(2);x(1)*x(2)-3*x(3)];
    y = C*x+dt*R^2*randn(n,1);

    Phi = [1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*x(2) y(1)*y(3) y(2)*y(3) sin(y(1)) sin(y(2)) sin(y(3)) cos(y(1)) cos(y(2)) cos(y(3)) zeros(32,1)';
           zeros(16,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*x(2) y(1)*y(3) y(2)*y(3) sin(y(1)) sin(y(2)) sin(y(3)) cos(y(1)) cos(y(2)) cos(y(3)) zeros(16,1)';
           zeros(32,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*x(2) y(1)*y(3) y(2)*y(3) sin(y(1)) sin(y(2)) sin(y(3)) cos(y(1)) cos(y(2)) cos(y(3))];
    
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
plot3(yArray(1,:),yArray(2,:),yArray(3,:),'-b','LineWidth',6);
hold on;
plot3(xbarArray(1,:),xbarArray(2,:),xbarArray(3,:),':r','LineWidth',6)
legend('measured','estimated')
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
subplot(3,2,1)
plot(t,sigma*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(2,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true','estimated')
grid on;
grid minor;
ylim([0 20])
ylabel('\sigma')
subplot(3,2,2)
plot(t,-sigma*ones(1,length(t))-thetabarArray(2,:)'/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \sigma')
subplot(3,2,3)
plot(t,rho*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(18,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([25 30])
ylabel('\rho')
subplot(3,2,4)
plot(t,rho*ones(1,length(t))-thetabarArray(18,:)'/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \rho')
subplot(3,2,5)
plot(t,beta*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(36,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylim([0 6])
ylabel('\beta')
xlabel('t (s)')
subplot(3,2,6)
plot(t,-beta*ones(1,length(t))-thetabarArray(36,:)'/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('error \beta')
xlabel('t (s)')

Coeff = round([(1/dt)*thetabar(1:(r/n),end)'; (1/dt)*thetabar((r/n)+1:2*(r/n),end)'; (1/dt)*thetabar(2*(r/n)+1:r,end)'])