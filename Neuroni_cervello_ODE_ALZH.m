% Michiel Bertsch, Bruno Franchi, Luca Meacci, Mario Primicerio, and Maria Carla Tesi
% The amyloid cascade hypothesis and Alzheimer's disease: a mathematical model
% European Journal of Applied Mathematics, 2020
% ---
% SCRIPT Figure 5

Tmax=40;
dim=20;

%Variables
m=10^(-2);
kappa=10^(-4);
lambda0=2;
tSpan=[0 180];
x0=zeros(dim,dim);
y0=20*ones(dim,dim);
z0=zeros(dim,dim);
beta=15;
theta=10^(-3);
ystar=22;
lambda0=2;
sigma=0.05;

a=zeros(dim,dim);
a(6,6)=0.02;
lambda=lambda0*ones(dim,dim);
yall=zeros(dim,dim);


for i=1:Tmax
    for j=1:dim
        for k=1:dim
            yzero=[x0(j,k) y0(j,k) z0(j,k)];
options=odeset('RelTol',1E-11,'AbsTol',1E-11);
[t,y]=ode45(@(t,y)Sisdif_alzhm(t,y,m,lambda(j,k),kappa),tSpan,yzero,options);


yall(j,k)=y(end,2);

x0(j,k)=y(end,1);
y0(j,k)=y(end,2);
z0(j,k)=y(end,3);

a(j,k)=a(j,k)+theta*subplus(y(end,2)-ystar);

        end
    end
%Update
ag=[zeros(1,dim);a;zeros(1,dim)];
ag=[zeros(dim+2,1) ag zeros(dim+2,1)];    
    
for jj=1:dim
    for kk=1:dim
    Al=ag(jj+1,kk+1).*ones(3,3);
    Ab=ag(jj:jj+2,kk:kk+2);
    D=Ab-Al;
    D=(abs(D)+D)/2;
    a(jj,kk)=a(jj,kk)+sigma*sum(sum(D));
    lambda(jj,kk) = lambda0*(1-a(jj,kk))*(1+beta*a(jj,kk));
    end
end

t=0.5*i;


imagesc(a);
colormap(flipud(jet))
title(['a at year ' num2str(t) ' ' ]) 
colorbar;
caxis([0 1]);
saveas(gcf, ['SRS_time' ,num2str(i), '.png']);
%pause(0.1)

end

%%%%



