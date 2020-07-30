% Michiel Bertsch, Bruno Franchi, Luca Meacci, Mario Primicerio, and Maria Carla Tesi
% The amyloid cascade hypothesis and Alzheimer's disease: a mathematical model
% European Journal of Applied Mathematics, 2020
% ---
% SCRIPT Figure 7

lambda=2;
M=10^(-2);
kappa=10^(-4);
kappastar=kappa/20;
L=lambda/kappa;
a=M/kappa;
b=kappastar/M;
Xp=[];
Xp2=[];
Yp=1:30;
for i=1:30
    Y=i;
X=1./(2+2*b*Y)*(-(Y+a+1/2*b*Y^2)+sqrt((Y+a+1/2*b*Y^2)^2+4*L*(1+b*Y)));
Xp=[Xp X];
X2=Y+b*Y^2+sqrt((Y+b*Y^2)^2+2*Y^2+b*Y^3+2*a*Y);
Xp2=[Xp2 X2];
end
plot(Xp,'g','LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('Y')
hold on
plot(Xp2,'m','LineWidth',1.5)
set(gca,'FontSize',12)
xlabel('Y')
