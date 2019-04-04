function T=tridiag(alfa,beit,gama)
%------alfa peit gama are the same demontion vectors ;beit(1)=0 gama(1)=0]
%--------²úÉú3¶Ô½ÇÕó
%----A=diag(alfa)+diag(beit,1)+diag(gama,-1);
[m,n]=size(beit);
for i=1:n-1
    T(i,i+1)=gama(i+1);
end
for i=1:n
    T(i,i)=alfa(i);
end
for i=2:n
    T(i,i-1)=beit(i);
end