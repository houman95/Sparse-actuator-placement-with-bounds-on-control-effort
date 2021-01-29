function [delta,Wd] = alg_4(W,Wn,r,a0p,a0plim)
n=length(Wn);
Wtild_term = 2*max(eig(Wn));
Wntild = Wn/Wtild_term;
a0 = 1e-12;
c = 1e-6;
Etild = 30000;
e = min (0.5,exp(-Etild));
[delta_c,Wdc,~] = alg_3(W,Wn,Etild,c,a0);

[~,R,~]=qr(Wntild);
l = -real(sum(log(nonzeros(diag(R)))));

[~,R,~]=qr(Wdc);
u = -real(sum(log(nonzeros(diag(R)))));
    

u = 6000; %high upper bound for small selections of r

Etild = 0.5*(l+u);
e = min (0.5,exp(-Etild));
flag = 0;
while flag ~=1
    while u-l > a0p
        [delta,Wd] = alg_3(W,Wn,Etild,c,a0);
        if sum(delta) > r
            l = Etild;  Etild = (l+u)/2;
        else
            u = Etild; Etild = (l+u)/2;
        end
        e = min (0.5,exp(-Etild));
%         disp('Bound length=')
%         disp(u-l)
    end
    if sum(delta)~=r
        a0p=a0p/2;
        if a0p<a0plim
            break
        end
    else
        flag=1;
    end
end
end

