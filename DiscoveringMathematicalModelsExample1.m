%% Research code by Agus Hasan

clear;
clc;

%% time horizon
tf  = 5;        % time horizon
dt  = 0.001;    % time step
t   = dt:dt:tf; % time array

%% number of variables and coefficients
n = 2;          % number of state
r = 10;         % number of parameters 

%% system description
A = eye(n);
C = eye(n);

%% noise
R = 0.625;      % standard deviation of the noise

%% state initialization
x        = [10;0];
xbar     = x;
y        = x;
thetabar = zeros(r,1);
 
%% true parameters
m = 1;
k  = 10;
b  = 0.5;

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
Px      = 0.1*eye(n);
Pt      = 0.1*eye(r);
Gamma   = 1*zeros(n,r);

%% simulation
for i=1:(tf/dt)

    xArray         = [xArray x];
    xbarArray      = [xbarArray xbar];
    yArray         = [yArray y];
    thetabarArray  = [thetabarArray thetabar]; 
    
    x = A*x+dt*[0 1; -k/m -b/m]*x+dt*[0 1/m]';
    y = C*x+R*randn(n,1);

    Phi = [y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) zeros(5,1)';
           zeros(5,1)' y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2)];
    
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
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('x [m]','FontSize',72)
subplot(2,1,2)
plot(t,yArray(2,:),'-b','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('v [m/s]','FontSize',72)
xlabel('t (s)')

figure(2)
subplot(2,1,1)
plot(t,yArray(1,:),'-b','LineWidth',10);
hold on;
plot(t,xbarArray(1,:),':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('measured','estimated')
grid on;
grid minor;
ylabel('x [m]','FontSize',72)
subplot(2,1,2)
plot(t,yArray(2,:),'-b','LineWidth',10);
hold on;
plot(t,xbarArray(2,:),':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
ylabel('v [m/s]','FontSize',72)
xlabel('t (s)')

figure(3)
subplot(2,2,1)
plot(t,(0)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('true','estimated')
xlim([0 tf])
ylim([-2000 500])
grid on;
grid minor;
ylabel('\theta_1')
axes('Position',[.25 .62 .2 .2])
box on
plot(t,(0)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:)/dt,':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-1 1])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
subplot(2,2,2)
plot(t,(1)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(2,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
xlim([0 tf])
ylim([-2000 500])
grid on;
grid minor;
ylabel('\theta_2')
axes('Position',[.7 .62 .2 .2])
box on
plot(t,(1)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(2,:)/dt,':g','LineWidth',10);
xlim([tf-1 tf])
ylim([0 2])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
subplot(2,2,3)
plot(t,(-k/m)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(6,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
xlim([0 tf])
ylim([-2000 500])
grid on;
grid minor;
ylabel('\theta_6')
xlabel('t (s)')
axes('Position',[.25 .15 .2 .2])
box on
plot(t,(-k/m)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(6,:)/dt,':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-11 -9])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
subplot(2,2,4)
plot(t,(-b/m)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(7,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',36)
ylim([-2000 500])
xlim([0 tf])
grid on;
grid minor;
ylabel('\theta_7')
xlabel('t (s)')
axes('Position',[.7 .15 .2 .2])
box on
plot(t,(-b/m)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(7,:)/dt,':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-1 0])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;


Coeff = round([(1/dt)*thetabar(1:(r/n),end)'; (1/dt)*thetabar((r/n)+1:r,end)'],21)
