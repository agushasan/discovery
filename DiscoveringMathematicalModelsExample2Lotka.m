%% Research code by Agus Hasan
%% Paper: "WyNDA: A Method to Discover Mathematical Models of Dynamical Systems from Data" submitted for MethodsX.

clear;
clc;

%% time horizon
tf  = 20;           % time horizon
dt  = 0.01;         % time step
t   = dt:dt:tf;     % time array

%% number of variables and coefficients
n = 2;              % number of measured state
r = 14;             % number of estimated parameters 

%% noise
R = 0;              % noise covariance

%% state initialization
x        = [8;8];           % actual state initialization
xbar     = x;               % estimated state initialization
y        = x;               % measurement initialization
thetabar = zeros(r,1);      % estimated parameter initialization
 
%% true parameters
alpha = 1;
beta = 0.2;
delta = 0.1;
gamma = 0.2;

%% for plotting
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
    
    x = x+dt*[alpha*x(1)-beta*x(1)*x(2);delta*x(1)*x(2)-gamma*x(2)];
    y = x+dt*R^2*randn(n,1);

    Phi = [y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3 zeros(7,1)';
           zeros(7,1)' y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3];
    
    % Estimation using adaptive observer
    Kx = Px*inv(Px+Rx);
    Kt = Pt*Gamma'*inv(Gamma*Pt*Gamma'+Rt);
    Gamma = (eye(n)-Kx)*Gamma;

    xbar = xbar+(Kx+Gamma*Kt)*(y-xbar);
    thetabar = thetabar-Kt*(y-xbar);

    xbar = xbar+Phi*thetabar;

    thetabar = thetabar;
    Px = (1/lambdav)*eye(n)*(eye(n)-Kx)*Px*eye(n);
    Pt = (1/lambdat)*(eye(r)-Kt*Gamma)*Pt;
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
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlabel('t (s)')

figure(2)
subplot(4,2,1)
plot(t,alpha*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
legend('true','estimated')
grid on;
grid minor;
ylim([alpha-0.2 alpha+0.2])
ylabel('\alpha','FontSize',36)
subplot(4,2,2)
plot(t,alpha*ones(1,length(t))-thetabarArray(1,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \alpha','FontSize',36)
subplot(4,2,3)
plot(t,beta*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(5,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([0 0.8])
ylabel('\beta','FontSize',36)
subplot(4,2,4)
plot(t,beta*ones(1,length(t))+thetabarArray(5,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \beta','FontSize',36)
subplot(4,2,5)
plot(t,delta*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(12,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([0 0.2])
ylabel('\delta','FontSize',36)
subplot(4,2,6)
plot(t,delta*ones(1,length(t))-thetabarArray(12,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \delta','FontSize',36)
subplot(4,2,7)
plot(t,gamma*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(9,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([0 0.8])
ylabel('\gamma','FontSize',36)
xlabel('t (s)')
subplot(4,2,8)
plot(t,gamma*ones(1,length(t))+thetabarArray(9,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \gamma','FontSize',36)
xlabel('t (s)')
