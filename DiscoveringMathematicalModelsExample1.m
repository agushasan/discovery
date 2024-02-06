%% Research code by Agus Hasan
%% Paper: "WyNDA: A Method to Discover Mathematical Models of Dynamical Systems from Noisy Data" submitted for MethodsX.

clear;
clc;

load DATAMSD.mat;   % load measurement data from the actual system (mass-spring-damper system)

%% simulation horizon
tf  = 5;            % time horizon
dt  = 0.001;        % time step
t   = dt:dt:tf;     % time array

%% number of variables and coefficients
n = 2;              % number of measured state
r = 10;             % number of estimated parameters 

%% state initialization
x        = [DATAMSD(1,1);DATAMSD(2,1)];         % actual state initialization
xbar     = x;                                   % estimated state initialization
thetabar = zeros(r,1);                          % estimated parameter initialization

%% for plotting
xArray          = [];   % for plotting the actual state
xbarArray       = [];   % for plotting the estimated state
thetabarArray   = [];   % for plotting the estimated parameter

%% Initialization for estimator

lambdav = 0.995;        % define parameters for the adaptive observer algorithm
lambdat = 0.999;
Rx      = 1*eye(n);
Rt      = 1*eye(n);
Px      = 0.1*eye(n);
Pt      = 0.1*eye(r);
Gamma   = 1*zeros(n,r);

%% simulation
for i=1:(tf/dt)

    xArray         = [xArray x];                % data for plotting the actual state
    xbarArray      = [xbarArray xbar];          % data for plotting the estimated state
    thetabarArray  = [thetabarArray thetabar];  % data for plotting the estimated parameter
    
    y(1) = DATAMSD(1,i);                        % load the measurement data
    y(2) = DATAMSD(2,i);
    y = [y(1) y(2)]';

    Phi = [y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2) zeros(5,1)';       % define the approximation function
           zeros(5,1)' y(1) y(2) y(1)^2 y(2)^2 y(1)*y(2)];
    
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

%% Plotting the results

figure(1)
subplot(2,1,1)
plot(t,DATAMSD(1,:),'-b','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('x [m]','FontSize',24)
subplot(2,1,2)
plot(t,DATAMSD(2,:),'-b','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('v [m/s]','FontSize',24)
xlabel('t (s)')

figure(2)
subplot(2,1,1)
plot(t,DATAMSD(1,:),'-b','LineWidth',10);
hold on;
plot(t,xbarArray(1,:),':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
legend('measured','estimated')
grid on;
grid minor;
ylabel('x [m]','FontSize',24)
subplot(2,1,2)
plot(t,DATAMSD(2,:),'-b','LineWidth',10);
hold on;
plot(t,xbarArray(2,:),':r','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
ylabel('v [m/s]','FontSize',24)
xlabel('t (s)')

figure(3)
subplot(2,2,1)
plot(t,(0)*ones(1,length(t))*dt,'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:),':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
legend('true','estimated')
xlim([0 tf])
ylim([-2 0.5])
grid on;
grid minor;
ylabel('\theta_1')
axes('Position',[.35 .66 .1 .1])
box on
plot(t,(0)*ones(1,length(t))*dt,'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(1,:),':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-1 1])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
subplot(2,2,2)
plot(t,(1)*ones(1,length(t))*dt,'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(2,:),':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
xlim([0 tf])
ylim([-2 0.5])
grid on;
grid minor;
ylabel('\theta_2')
axes('Position',[.79 .66 .1 .1])
box on
plot(t,(1)*ones(1,length(t))*dt,'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(2,:),':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-0.005 0.005])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
subplot(2,2,3)
plot(t,(-84*dt)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(6,:),':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
xlim([0 tf])
ylim([-2 0.5])
grid on;
grid minor;
ylabel('\theta_6')
xlabel('t (s)')
axes('Position',[.35 .2 .1 .1])
box on
plot(t,(-84*dt)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(6,:),':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-0.1 -0.06])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;
subplot(2,2,4)
plot(t,(-0.9*dt)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(7,:),':g','LineWidth',10);
set(gca,'color','white','LineWidth',3,'FontSize',24)
ylim([-2 0.5])
xlim([0 tf])
grid on;
grid minor;
ylabel('\theta_7')
xlabel('t (s)')
axes('Position',[.79 .2 .1 .1])
box on
plot(t,(-0.9*dt)*ones(1,length(t)),'-k','LineWidth',10);
hold on;
plot(t,thetabarArray(7,:),':g','LineWidth',10);
xlim([tf-1 tf])
ylim([-0.001 -0.0006])
set(gca,'color','white','LineWidth',3,'FontSize',24)
grid on;
grid minor;

Coeff = round([thetabar(1:(r/n),end)'; thetabar((r/n)+1:r,end)'],3)
