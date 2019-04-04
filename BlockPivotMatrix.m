function [A]=BlockPivotMatrix(m)%[A]=BlockPivotMatrix(m,mu,eta,zeta)
%-----produce a block three diagonal matrix
A=tridiag(2*ones(1,m^2),[zeros(1,m^2)],[0 (-1)*ones(1,m^2-1)]);
i=1;
for i=1:m^2-m
    A(i,i+m)=-1;
    if i<m
   A(m*i,m*i+1)=0;
end
end
n=m^2; 
for i=1:n
    x(i)=(1/2)*(1+(-1)^(i-1));
end
x=x';
zo=2*ones(n,1)-x;
A=A+A';
%A=A+mu*eye(m^2,m^2)+eta*tridiag(zeros(1,m^2),[zeros(1,m^2)],[0 ones(1,m^2-1)])+zeta*diag(zo);
 
        