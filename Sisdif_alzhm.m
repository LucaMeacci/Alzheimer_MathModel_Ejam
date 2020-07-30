function dy=Sisdif_alzhm(t,y,m,lambda,k)
kstar=k/20;
dy= zeros(3,1); % vettore colonna
dy(1)= -k*y(1)^2-k*y(1)*y(2)-kstar*y(1)*y(3)-m*y(1)+lambda;
dy(2)= 1/2*k*y(1)^2-k*y(1)*y(2)-k*y(2)^2-kstar*y(2)*y(3)-m*y(2);
dy(3)=1/2*k*y(2)^2+k*y(1)*y(2)-m*y(3);
end