N = 2^10; % space discretization
M = 2^12; % time discretization
K = 2000; % coefficients


x = linspace(0, 1, N);
t = linspace(0, 1, M);
%% 
y = zeros(M, N);

coef = 1.0 ./ (1:K).^2 .* (-1).^(1:K);

tic;
for time_step = 1:M
    for i = 1:K
        y(time_step, :) = y(time_step, :) + coef(i) * exp(-pi^2 * i^2 * t(time_step)) * sin( i * pi * x );
    end
end
toc;
%%

tic;
[U_, S_, V_] = svd(y);
toc;
%%
V = V_(:,1:20);
tau_ = 1e-3:1e-3:0.1;

int_err = zeros(size(tau_, 2), 1);

for t=1:size(tau_, 2)
tau = tau_(t);
S = 50;
% fprintf('sample size %d\n', S);
H = 1000; % number of experiments
sample_id = 1:ceil(N/S):N;
err = 0;

% K = ceil( sqrt(-log(10^(-12)) / pi^2 / tau ) );
K = 2000;
% tic;
for p = 1:H
    f_coef = 2 * rand(1, K) - 1; % all between -1 and 1.
    f_coef = f_coef / norm(f_coef);

    z_true = zeros(1, N);
    z_samp = zeros(1, S);
    
    for i = 1:K
        z_true = z_true + f_coef(i) * exp(-pi^2 * i^2 * tau) * sin(i * pi * x);
    end
    
    z_sample = z_true(sample_id);
    coeff = V(sample_id, :) \ z_sample';
    z_rec = V * coeff;
    
%     fprintf('p=%d, error=%e\n',...
%         p,  norm(z_rec - z_true')/norm(z_true) );
    
    err =  err + norm(z_rec - z_true')/norm(z_true);
end
% toc;
fprintf('tau=%f, average error=%e\n', tau, err / H);
int_err(t) = err/H;
end



