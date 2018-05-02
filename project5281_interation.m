%syms x;
%initial polynomial p
% rp=[-1.8337 1.7696];
rp=[-1.8337];
p=poly(rp);
P=poly2sym(p);
%initial polynomial f
f=zeros(1,6);
f(1)=1;
% rf=[3.42 3.46 3.5 3.54 3.58];
% f=poly(rf);
% F=poly2sym(f);
%derivation of c function
S=[-1 -0.6 -0.2 0.2 0.6 1];
for num=1:5
    
    %solve the linear equations
    M=[S(1)^4 S(1)^3 S(1)^2 S(1) 1 polyval(p,S(1)),
       S(2)^4 S(2)^3 S(2)^2 S(2) 1 -polyval(p,S(2)),   
       S(3)^4 S(3)^3 S(3)^2 S(3) 1 polyval(p,S(3)),
       S(4)^4 S(4)^3 S(4)^2 S(4) 1 -polyval(p,S(4)),
       S(5)^4 S(5)^3 S(5)^2 S(5) 1 polyval(p,S(5)),
       S(6)^4 S(6)^3 S(6)^2 S(6) 1 -polyval(p,S(6))];
    V=-[S(1)^5 S(2)^5 S(3)^5 S(4)^5 S(5)^5 S(6)^5];
    X=inv(M)*V';
    for k=1:5
        f(k+1)=X(k);
    end
    k=0;
    F=poly2sym(f);
    
    [df,dp]=polyder(f,p);
    z=roots(df);
    %choose inner passband extreme point
    for k=1:length(z)
        if z(k)>-1 && z(k)<1
           z(k)=z(k);
        else 
           z(k)=-2;
        end
    end
    k=0;
    %order those point
    z=sort(z);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%num of Tz
%     for k=1:length(z)-2
%         S(k+1)=z(k+2);
%     end
    for k=1:length(z)
        if z(k)==-2
           S(k)=S(k);
        else S(k)=z(k);
        end
    end
    k=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%num of Tz    
            
end
%after interation
C=F/(P*X(6));
    fplot(C,[-1,1],'r');
    line([-1,1],[0,0],'linestyle','--');