% Michiel Bertsch, Bruno Franchi, Luca Meacci, Mario Primicerio, and Maria Carla Tesi
% The amyloid cascade hypothesis and Alzheimer's disease: a mathematical model
% European Journal of Applied Mathematics, 2020
% ---
% SCRIPT Figure 1

%Variables
m=10^(-2);
k=10^(-4);
lambda=2;
tSpan=[0 500];
y0 =[0 0 0];
options=odeset('RelTol',1E-11,'AbsTol',1E-11);
[t,y]=ode45(@(t,y)Sisdif_alzhm(t,y,m,lambda,k),tSpan,y0,options);

x=y(:,1);
w=y(:,2);
z=y(:,3);

%Plotting
figure(1)
plot(t,y,'LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('time (days)')

%%%%

w(end)

