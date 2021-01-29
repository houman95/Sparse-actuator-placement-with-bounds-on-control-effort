clear;
n=5;
A = [-1 0 0 0 0;
     1 -1 0 0 0;
     0 1 -1 0 0;
     0 0 1 -1 0;
     0 0 0 1 -1];
W = zeros(n,n,n);
for i = 1:n
    B = zeros(n);
    B(i,i)=1;
    sys = ss(A,B,zeros(1,n),0);
    W(:,:,i) = gram(sys,'c');
end
B = eye(5);
sys = ss(A,B,zeros(1,n),0);
Wn = zeros(n);
for i = 1:n
    Wn = Wn + W(:,:,i);
end
Wtild_term = 2*max(eig(Wn));
Etild_term = n*log(2*max(eig(Wn)));
Wntild = Wn/Wtild_term;
%% Step 2
    % E >= log(det(inv(Wn)))
    % Etild >= log(det(inv(Wn))) + Etild_term
a0 = 1e-6;
c = 1e-6;
Etildmin = round(log(det(inv(Wn))) + Etild_term);
Etild_v = linspace(Etildmin,20+Etildmin,11);
delta_v = zeros(n,10);
e_v = zeros(n,1);
for i=1:11
    [delta_v(:,i),~,e_v(i)] = alg_3(W,Wn,Etild_v(i),c,a0);
end
%% Step 3
Wd_logdet = zeros(n,1);
delta_v = zeros(n,n);
for r = 1:n
    [delta_r,Wd_r] = alg_4(W,Wn,r,10,1e-12);
    if sum(delta_r)>r
        disp('error in algorithm 4')
%         break
    end
    delta_v(:,r) = delta_r; %saves every delta for every r
    [~,R,~]=qr(Wd_r);
    sumwd = real(sum(log(nonzeros(diag(R))))); %sumwd = log(det(W_delta_r))
    Wd_logdet (r) = sumwd;
    Wd_det (r) = det(Wd_r);
end
figure(1)
plot(1:n,Wd_logdet)
xlabel('r')
ylabel('Log(Det(W_\Delta))')
grid on
%% Step 4
%% Initialization
load('GD99_c.mat');
load('alpha.mat');
load('beta.mat');
S = Problem.A;
n = 105;
% alpha = rand(n,1);
% beta = rand(n,1);
for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = -alpha(i) + beta(i)*S(i,j);
        else
            A(i,j) = beta(i)*S(i,j);
        end
    end
end
W = zeros(n,n,n);
for i = 1:n
    B = zeros(n);
    B(i,i)=1;
    W(:,:,i) = CtrGram(A,B);
    sys = ss(A,B,zeros(1,n),0);
%     W(:,:,i) = gram(sys,'c');
end
% W = round(W,4);
Wn = zeros(n);
for i = 1:n
    Wn = Wn + W(:,:,i);
end
Wtild_term = 2*max(eig(Wn));
Etild_term = n*log(2*max(eig(Wn)));
Wntild = Wn/Wtild_term;
%% Step 2 for Step 4
a0 = 1e-6;
c = 1e-6;
E = 10;
Etild = E + Etild_term;
[delta,Wd,e] = alg_3(W,Wn,Etild,c,a0);
deltaset =[];
for i = 1:n
    if delta(i)==1
        deltaset=[deltaset i]; % set of actuators
    end
end
%% Step 3 for Step 4
Wd_logdet = zeros(n,1);
delta_v = zeros(n,n);
for r = 1:n
    r
    [delta_r,Wd_r] = alg_4(W,Wn,r,100,1);
    if sum(delta_r)>r
        disp('error in algorithm 4')
%         break
    end
    delta_v(:,r) = delta_r; %saves every delta for every r
    [~,R,~]=qr(Wd_r);
    sumwd = real(sum(log(nonzeros(diag(R))))); %sumwd = log(det(W_delta_r))
    Wd_logdet (r) = sumwd;
    Wd_det (r) = det(Wd_r);
end
%%
load ('delta_v.mat'); % We already saved delta_v for our system so we can comment the procedure above
Wd_logdet=zeros(n,1);
Wd_det = zeros(n,1);
for k=1:n
    Wd = zeros(n);
    for i = 1:n
        Wd = Wd+delta_v(i,k)*W(:,:,i); %forming W_Delta
    end
     [~,R,~]=qr(Wd);
    sumwd = real(sum(log(nonzeros(diag(R)))));
    Wd_logdet(k)=sumwd;
    Wd_det(k) = det(Wd);
end
figure(2)
plot(1:n,Wd_logdet)
xlabel('r')
ylabel('Log(Det(W_\Delta))')
grid on
%% Investigating Optimality
Wd5 = Wd_logdet(5);
Wd_logdet_set = zeros(1000,1);
for i = 1:1000
    randset = randperm(105,5); % choose 5 random unique integers from 1 to 105
    Wd = W(:,:,randset(1))+W(:,:,randset(2))+W(:,:,randset(3))+...
        W(:,:,randset(4))+W(:,:,randset(5)); %forming W_Delta
    [~,R,~]=qr(Wd);
    sumwd = real(sum(log(nonzeros(diag(R))))); % sumwd = log(det(W_delta))
    Wd_logdet_set(i) = sumwd;
end
figure(3)
hist(Wd_logdet_set,100)
xlabel('log(det(W_\Delta))')
ylabel('count')
axis([-6000 -2000 0 32])
hold on 
stem(Wd5,32,'r')
legend('histogram for random \Delta','log(det(W_\Delta)) obtained by our algorithm')