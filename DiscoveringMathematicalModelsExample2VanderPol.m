%% Research code by Agus Hasan
%% Paper: "WyNDA: A Method to Discover Mathematical Models of Dynamical Systems from Data" submitted for MethodsX.

clear;
clc;

%% time horizon
tf  = 5;            % time horizon
dt  = 0.001;        % time step
t   = dt:dt:tf;     % time array

%% number of variables and coefficients
n = 2;              % number of measured state
r = 18;             % number of estimated parameters 

%% noise
R = 0;              % noise covariance

%% state initialization
x        = [1;1];           % actual state initialization
xbar     = x;               % estimated state initialization
y        = x;               % measurement initialization
thetabar = zeros(r,1);      % estimated parameter initialization
 
%% true parameters
mu = 1.2;

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
    
    x = x+dt*[x(2);mu*(1-x(1)^2)*x(2)-x(1)];
    y = x+dt*R^2*randn(n,1);

    Phi = [y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3 y(1)^2*y(2) y(1)*y(2)^2 zeros(9,1)';
           zeros(9,1)' y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) y(1)^3 y(2)^3 y(1)^2*y(2) y(1)*y(2)^2];
    
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
%legend('true # of prey','estimated # of prey','true # of predator','estimated # of predator')
grid on;
grid minor;
%ylabel('# of population')
xlabel('t (s)')

figure(2)
subplot(2,1,1)
plot(t,mu*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(11,:)/dt,':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
legend('true','estimated')
grid on;
grid minor;
ylim([mu-1 mu+1])
ylabel('\mu','FontSize',36)
subplot(2,1,2)
plot(t,mu*ones(1,length(t))-thetabarArray(11,:)/dt,':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('error \mu','FontSize',36)
xlabel('t (s)')

figure(3)
plot(yArray(1,:),yArray(2,:),'-','LineWidth',16);
hold on;
plot(xbarArray(1,:),xbarArray(2,:),':','LineWidth',16)
set(gca,'color','white','LineWidth',3,'FontSize',36)
legend('measured','estimated')
grid on;
grid minor;
ylabel('q')
xlabel('p')
