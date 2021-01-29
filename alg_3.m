function [delta,Wd,e] = alg_3(W,Wn,Etild,c,a0)
n = length(Wn);
Wtild_term = 2*max(eig(Wn));
a = a0;
flag = 0;
l=0;
u=min (0.5,exp(-Etild));
e = 0.5*(l+u);
if e==0
    flag=1;
end
[delta,~] = alg_2(W,Wn,Etild,e);
Wd = zeros(n);
for i = 1:n
   Wd = Wd+delta(i)*W(:,:,i);
end
Wdtild = Wd/Wtild_term;
while flag~=1
   while u-l > a
       [delta,~] = alg_2(W,Wn,Etild,e);
       Wd = zeros(n);
       for i = 1:n
           Wd = Wd+delta(i)*W(:,:,i);
       end
       Wdtild = Wd/Wtild_term;
       if -log(det(Wdtild))+log(det(Wdtild+e*eye(n))) > c*Etild
           u=e;
       else
           l=e;
       end
       e = 0.5*(l+u);
   end
   if -log(det(Wdtild))+log(det(Wdtild+e*eye(n))) > c*Etild
       u=e; e=0.5*(l+u);
   end
   [delta,~] = alg_2(W,Wn,Etild,e);
   if -log(det(Wdtild))+log(det(Wdtild+e*eye(n))) <= c*Etild
       flag = 1;
   else
       a = a/2;
   end
end
end

