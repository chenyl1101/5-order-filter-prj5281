Mr=M;%reserve initial M
%rotation matrix R
R=eye(7,7);
%annihilate element m,n 0-->S 1-5 7-->L
%inner(m,n)-->M(m+1,n+1)
% m=0;n=5;
for m=0:3
    n=5;
    while n-m>=2
%pivotinner(n-1,n)-->pivotM(n,n+1)=pivotM(i,j)
%calculate theta:Mkl/Mmn-->k=m+1,l=n+1,m=m+1,n=n;
          theta_r=-atan(Mr(m+1,n+1)/Mr(m+1,n));
          cr=cos(theta_r);
          sr=sin(theta_r);
          R(n,n)=cr;
          R(n+1,n+1)=cr;
          R(n+1,n)=sr;
          R(n,n+1)=-sr;
          Mr=R*Mr*R';
          R=eye(7,7);
          n=n-1;
    end
end
%arrow topology obtained
%creat one tri-section to contribute TZ/rp
%1st rotarion: pivotinner(4,5)
m=4;n=5;w=rp;
theta_r=atan(Mr(m+1,n+1)/(w+Mr(n+1,n+1)));
cr=cos(theta_r);
sr=sin(theta_r);
R(n,n)=cr;
R(n+1,n+1)=cr;
R(n+1,n)=sr;
R(n,n+1)=-sr;
Mr=R*Mr*R';
R=eye(7,7);
%2nd rotation: pivotinner(3,4)
m=3;n=4;w=rp;
theta_r=atan(Mr(m+1,n+1)/(w+Mr(n+1,n+1)));
cr=cos(theta_r);
sr=sin(theta_r);
R(n,n)=cr;
R(n+1,n+1)=cr;
R(n+1,n)=sr;
R(n,n+1)=-sr;
Mr=R*Mr*R';
R=eye(7,7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set main coupling be positive
k=0;
for k=1:6
    if Mr(k,k+1)<0
       Mr(:,k+1)=-Mr(:,k+1);
       Mr(k+1,:)=-Mr(k+1,:);
    end
end
%verify coulping matrix M
Gm=zeros(7,7);
Gm(1,1)=1;
Gm(7,7)=1;
Cm=eye(7,7);
Cm(1,1)=0;
Cm(7,7)=0;
syms x;
Zm=Gm+i*x.*Cm+i.*Mr;% MorMr
Zm=inv(Zm);
s21=2*Zm(7,1);
s11=-1+2*Zm(1,1);
mags11=10*log(abs(s11));
mags21=10*log(abs(s21));
fplot(mags11,[-5,5],'b');hold on;
fplot(mags21,[-5,5],'r');
