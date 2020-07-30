% Michiel Bertsch, Bruno Franchi, Luca Meacci, Mario Primicerio, and Maria Carla Tesi
% The amyloid cascade hypothesis and Alzheimer's disease: a mathematical model
% European Journal of Applied Mathematics, 2020
% ---
% SCRIPT Figures 2, 3, and 4

Tmax=40;

%Variables
m=10^(-2);
k=10^(-4);
lambda=2;
tSpan=[0 180];
y0 =[0 20 0];
beta=15;
theta=10^(-3);
ystar=22;
lambda0=2;
lambda=lambda0;
a0=0.02;
a=a0;
as=a;

ys=y0(1,2);
lambdas=lambda;


for i=1:Tmax
options=odeset('RelTol',1E-11,'AbsTol',1E-11);
[t,y]=ode45(@(t,y)Sisdif_alzhm(t,y,m,lambda,k),tSpan,y0,options);

wend=y(end,2);
ys=[ys wend];
y0=y(end,:);

a=a(1,end)+theta*subplus(wend-ystar);
lambda = lambda0*(1-a)*(1+beta*a);
lambdas=[lambdas lambda];
as = [as a];

end

anni=[0:Tmax]/2;

%Plotting
figure(1)
plot(anni(2:end),as(2:end),'o black','LineWidth',1.5)
set(gca,'FontSize',12)
axis([0 Tmax/2 0 1])
xlabel('time (years)')
ylabel('a')

figure(2)
plot(anni(2:end),lambdas(2:end),'o','LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('time (years)')

figure(3)
plot(anni(3:end),ys(3:end),'r','LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('time (years)')
ylabel('y')

%%%%




