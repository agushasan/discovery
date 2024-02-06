%% Research code by Agus Hasan
%% Paper: "WyNDA: A Method to Discover Mathematical Models of Dynamical Systems from Data" submitted for MethodsX.

clear;
clc;

%% time horizon
tf  = 20;           % time horizon
dt  = 0.001;        % time step
t   = dt:dt:tf;     % time array

%% number of variables and coefficients
n = 3;              % number of measured state
r = 16;             % number of estimated parameters 

%% noise
R = 0;              % noise covariance

%% state initialization
x        = [-8;5;10];       % actual state initialization
xbar     = x;               % estimated state initialization
y        = x;               % measurement initialization
thetabar = zeros(r,1);      % estimated parameter initialization
 
%% true parameters
a = 0.1;
b = 0.1;
c = 14;

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
    
    x = x+dt*[-x(2)-x(3);x(1)+a*x(2);b+x(3)*(x(1)-c)];
    y = x+dt*R^2*randn(n,1);

    Phi = [y(1) y(2) y(3) zeros(13,1)';
           zeros(3,1)' y(1) y(2) y(3) zeros(10,1)';
           zeros(6,1)' 1 y(1) y(2) y(3) y(1)^2 y(2)^2 y(3)^2 y(1)*y(2) y(1)*y(3) y(2)*y(3)];
    
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
set(gca,'color','white','LineWidth',3,'FontSize',36)
grid on;
grid minor;
xlabel('p')
ylabel('q')
zlabel('r')

figure(2)
subplot(3,2,1)
plot(t,a*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(5,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
legend('true','estimated')
grid on;
grid minor;
ylim([0 0.2])
ylabel('a','FontSize',36)
subplot(3,2,2)
plot(t,-a*ones(1,length(t))-thetabarArray(5,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error a','FontSize',36)
subplot(3,2,3)
plot(t,b*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(7,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([0 0.2])
ylabel('b','FontSize',36)
subplot(3,2,4)
plot(t,b*ones(1,length(t))-thetabarArray(7,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error b','FontSize',36)
subplot(3,2,5)
plot(t,c*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,-thetabarArray(10,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylim([10 20])
ylabel('c','FontSize',36)
xlabel('t (s)')
subplot(3,2,6)
plot(t,-c*ones(1,length(t))-thetabarArray(10,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error c','FontSize',36)
xlabel('t (s)')
