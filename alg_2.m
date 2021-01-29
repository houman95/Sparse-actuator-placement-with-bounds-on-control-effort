function [delta,Wd] = alg_2(W,Wn,Etild,e)
n = length(Wn);
Wtild_term = 2*max(eig(Wn));
delta = zeros(n,1);
Nu = ones(n,1);
Wd = zeros(n);
Wdtild = Wd/Wtild_term;
a=0;
done = 0;
sumwd = 1e10; %Initialization with a high value
while (sumwd > Etild) && (done == 0)
    if a~=0
        delta(a)=1;
    end
    Nu_delta = Nu - delta;
    Wda = zeros(n,n,n);
    Wd = zeros(n);
    for i = 1:n
        Wd = Wd+delta(i)*W(:,:,i); %forming W_Delta
    end
    Wdtild = Wd/Wtild_term;
    [~,R,~]=qr(Wdtild+e*eye(n));
    sumwd = -real(sum(log(nonzeros(diag(R))))); %sumwd = log(det(W_delta tild))
    if sumwd==0 
        sumwd = 1e10; %to enter the loop again
    end
    for i = 1:n
        Wda(:,:,i) = Wd + Nu_delta(i)*W(:,:,i);
        Wdatild = Wda(:,:,i)/Wtild_term;
    [~,R,~]=qr(Wdatild+e*eye(n));
    sumwda(i) = -real(sum(log(nonzeros(diag(R)))));% sumwd = log(det(W_delta+a tild))
        logset(i) =  sumwd - sumwda(i); %saving all the values to choose 'a' that maximizes it
    end
    [~,a]=max(logset);
    if sum(Nu_delta)==0 %if all the actuators are chosen but we're still in the loop
        done = 1;
    end
end
end