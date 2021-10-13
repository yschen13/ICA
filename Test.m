% This script aims to test the performance of different ICA algorithms


% Test of ACMNSym 
ACMNsym(C_highspeed,'mle_noncirc');


% Synthetic dataset
N = 5; T = 20000;
% circular sources
rng(1); s_r = normrnd(0,1,[1 T]);
rng(2); s_i = normrnd(0,1,[1 T]);
s_circ1 = s_r + sqrt(-1)*s_i;
rng(3); s_r = random('Uniform',0,1,1,T);
rng(4); s_i = random('Uniform',0,1,1,T);
s_circ2 = s_r + sqrt(-1)*s_i;
Pcov1 = s_circ1*transpose(s_circ1);
Cov1 = s_circ1*s_circ1';
Pcov2 = s_circ2*transpose(s_circ2);
Cov2 = s_circ2*s_circ2';
% non-circular sources
rng(5); s_r = normrnd(0,1,[1 T]);
rng(6); s_i = normrnd(0,1,[1 T]);
s_noncirc1 = 2*s_r + sqrt(-1)*s_i;
rng(7); s_r = random('Uniform',0,1,1,T);
rng(8); s_i = random('Uniform',0,1,1,T);
s_noncirc2 = 2*s_r + sqrt(-1)*s_i;
rng(9); s_r = random('Rayleigh',1,1,T); % all positive values
rng(10); s_i = random('Uniform',0,1,1,T);
s_noncirc3 = s_r + sqrt(-1)*s_i;

% mix of circular and noncircular sources
s = [s_circ1; s_circ2; s_noncirc1; s_noncirc2; s_noncirc3];
rng(11)
A = random('Uniform',0,1,N,N); % mixing matrix
x = A*s; % observed samples
[a,W,WOr,alphas] = ACMNsym(x,'mle_noncirc');
C = W*A;
% PI = 4: OK performance


% all circular sources
rng(1); tmp = random('Uniform',0,1,[2*N T]);
s = tmp(1:N,:) + sqrt(-1)*tmp(N+1:end,:);
x = A*s; % observed sample
[a,W,WOr,alphas] = ACMNsym(x,'mle_circ');
C = W*A;
% PI = 0.6: very low/GOOD

% all noncircular sources
rng(1); tmp = random('Uniform',0,1,[2*N T]);
s = 3*tmp(1:N,:) + sqrt(-1)*tmp(N+1:end,:);
x = A*s; % observed sample
[a,W,WOr,alphas] = ACMNsym(x,'mle_noncirc');
C = W*A;
% PI = 0.6: very low/GOOD

% performance index
% C = eye(N); C = C(randperm(N),:); % PI of permutation matrix is 0
PI_p1 = sum(sum(abs(C),2) ./ max(abs(C),[],2) - 1);
PI_p2 = sum(sum(abs(C),1) ./ max(abs(C),[],1) - 1);
PI = PI_p1 + PI_p2