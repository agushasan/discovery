%% Research code by Agus Hasan
%% Paper: "WyNDA: A Method to Discover Mathematical Models of Dynamical Systems from Data" submitted for MethodsX.

clear;
clc;

%% time horizon
tf  = 5;            % time horizon
dt  = 0.001;        % time step
t   = dt:dt:tf;     % time array

%% number of variables and coefficients
n = 3;              % number of measured state
r = 30;             % number of estimated parameters 

%% noise
R = 1;              % noise covariance

%% state initialization
x        = [-8;7;27];       % actual state initialization
xbar     = x;               % estimated state initialization
y        = x;               % measurement initialization
thetabar = zeros(r,1);      % estimated parameter initialization
 
%% true parameters
sigma   = 10;
rho     = 28;
beta    = 3;

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
    
    x = x+dt*[10*(x(2)-x(1));x(1)*(28-x(3))-x(2);x(1)*x(2)-3*x(3)];
    y = x+dt*R^2*randn(n,1);

    Phi = [1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3) zeros(20,1)';
           zeros(10,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3) zeros(10,1)';
           zeros(20,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3)];
    
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
plot3(yArray(1,:),yArray(2,:),yArray(3,:),'-b','LineWidth',16);
hold on;
plot3(xbarArray(1,:),xbarArray(2,:),xbarArray(3,:),':r','LineWidth',16)
legend('measured','estimated')
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
xlabel('p')
ylabel('q')
zlabel('r')

figure(2)
subplot(3,2,1)
plot(t,sigma*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(2,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
legend('true','estimated')
grid on;
grid minor;
ylim([0 20])
ylabel('\sigma','FontSize',36)
subplot(3,2,2)
plot(t,-sigma*ones(1,length(t))-thetabarArray(2,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \sigma','FontSize',36)
subplot(3,2,3)
plot(t,rho*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(12,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([25 30])
ylabel('\rho','FontSize',36)
subplot(3,2,4)
plot(t,rho*ones(1,length(t))-thetabarArray(12,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \rho','FontSize',36)
subplot(3,2,5)
plot(t,beta*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(24,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([0 6])
ylabel('\beta','FontSize',36)
xlabel('t (s)')
subplot(3,2,6)
plot(t,-beta*ones(1,length(t))-thetabarArray(24,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \beta','FontSize',36)
xlabel('t (s)')
