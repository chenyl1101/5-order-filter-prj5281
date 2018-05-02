rootf=i*roots(f);
rootp=i*roots(p);%s=jw
Fs=poly(rootf);
Ps=poly(rootp); 
Ps=i*Ps;%N-nfz is even
% Fs=real(Fs);
% Ps=real(Ps);
E=abs(polyval(p,[1]))/abs(polyval(f,[1]))/sqrt(10^2.2-1);
%compute Es
rE=roots(sym2poly(P/E-i*F));
for k=1:length(rE)
    if imag(rE(k))<0
       rE(k)=conj(rE(k));
    else
       rE(k)=rE(k);
    end
end
rootE=i*rE;%s=jw
Ew=poly(rE);
Es=poly(rootE);
%compute S11 S21 with w
s11=F/poly2sym(Ew);
s21=P/(E*poly2sym(Ew));
mags11=10*log(abs(s11));
mags21=10*log(abs(s21));
fplot(mags11,[-5,5],'r');hold on;
fplot(mags21,[-5,5],'b');
% S11=poly2sym(Fs)/poly2sym(Es);
% S21=poly2sym(Ps)/(E*poly2sym(Es));

%since roots of Fs are on imag axis Fs*=Fs N=5
F22S=-Fs;
%convert to Y matrix
yd=zeros(1,6);
for k=1:6
    if mod(k,2)==0
       yd(k)=i*imag(Es(k)+Fs(k));
    else
       yd(k)=real(Es(k)+Fs(k));
    end
end
k=0;
y11n=zeros(1,5);
y22n=zeros(1,5);
for k=1:5
    if mod(k,2)==0
       y11n(k)=i*imag(Es(k+1)-Fs(k+1));
       y22n(k)=i*imag(Es(k+1)+Fs(k+1));
    else
       y11n(k)=real(Es(k+1)-Fs(k+1));
       y22n(k)=real(Es(k+1)+Fs(k+1));
    end
end
k=0;
y21n=-Ps/E;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yd=yd/2;y11n=y11n/2;y22n=y22n/2;y21n=y21n/2;
%y11=y11n/yd;y12=y21=y21n/yd;y22=y22n/yd
[r11k,pk,k11]=residue(y11n,yd);
[r21k,pk,k21]=residue(y21n,yd);
lamdak=pk/i;
r22k=r11k;
%N+2 transversal matrix M
M=zeros(7,7);
%M(1,1)=k11=0;M(7,7)=k22=0;M17=M71=k21=0;
for k=2:6
    M(k,k)=-lamdak(k-1);
    M(k,7)=sqrt(r22k(k-1));
    M(7,k)=M(k,7);
    M(1,k)=r21k(k-1)/sqrt(r22k(k-1));
    M(k,1)=M(1,k);
    
end
k=0;M=real(M);


