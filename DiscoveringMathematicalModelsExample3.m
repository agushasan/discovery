%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 1;
dt  = 0.0001;
t   = dt:dt:tf;

%% number of variables and coefficients
n = 2;
r = 14;

%% system description
A = eye(n);
C = eye(n);

%% noise
RF = 1;

%% state initialization
x        = [1;1];
xbar     = x;%zeros(n,1);
y        = x;%zeros(n,1);
thetabar = zeros(r,1);
 
%% true parameters
J = 0.01;
b = 0.1;
K = 0.2;
R = 1;
L = 0.5;

%% initial control inputs
u     = 2;

%% for plotting
uArray          = [];
xArray          = [];
xbarArray       = [];
yArray          = [];
thetabarArray   = [];

%% Initialization for estimator

lambdav = 0.995;
lambdat = 0.999;
Rx      = 1*eye(n);
Rt      = 1*eye(n);
Px      = 1*eye(n);
Pt      = 1*eye(r);
Gamma   = 1*zeros(n,r);

%% simulation
for i=1:(tf/dt)
    
    %u = 4*sin(i*dt*40);

    uArray         = [uArray u];
    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 
    
    x = A*x+dt*[-b/J K/J; -K/L -R/L]*x+dt*[0 1/L]'*u;
    y = C*x+dt*RF^2*randn(n,1);

    Phi = [y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) u u^2 zeros(7,1)';
           zeros(7,1)' y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) u u^2];
    
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
subplot(2,1,1)
plot(t,yArray(1,:),'-b','LineWidth',10);
hold on;
plot(t,xbarArray(1,:),':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('measured','estimated')
grid on;
grid minor;
ylabel('\omega [rad/s]')
subplot(2,1,2)
plot(t,yArray(2,:),'-b','LineWidth',10);
hold on;
plot(t,xbarArray(2,:),':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('i [Amp]')
xlabel('t (s)')

figure(2)
subplot(3,2,[1 2])
plot(t,uArray,'-k','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('control input')
xlim([0 tf])
grid on;
grid minor;
ylabel('u','FontSize',72)
xlabel('t (s)')
subplot(3,2,3)
plot(t,(-b/J)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true','estimated')
xlim([0 tf])
grid on;
grid minor;
ylabel('\theta_1','FontSize',72)
subplot(3,2,4)
plot(t,(K/J)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(2,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
xlim([0 tf])
grid on;
grid minor;
ylabel('\theta_2','FontSize',72)
subplot(3,2,5)
plot(t,(-K/L)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(8,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
xlim([0 tf])
grid on;
grid minor;
xlabel('t (s)')
ylabel('\theta_8','FontSize',72)
subplot(3,2,6)
plot(t,(-R/L)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(9,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
xlim([0 tf])
grid on;
grid minor;
xlabel('t (s)')
ylabel('\theta_9','FontSize',72)


Coeff = round([(1/dt)*thetabar(1:(r/n),end)'; (1/dt)*thetabar((r/n)+1:r,end)'],1)
