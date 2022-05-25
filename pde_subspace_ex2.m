N = 2^10;
M = 2^12;
K = 2000;

x = linspace(0, 1, N);
t = linspace(0, 1, M);

%%
y1 = zeros(M, N);
y2 = zeros(M, N);

coef1 = 1.0 ./ (1:K).^2 .* (-1).^(1:K);
coef2 = 1.0 ./ (1:K).^2 .* (-1).^(1:K);

eigen = pi^2 .* (1:K).^2 + 1;

for time_step = 1:M
    for i = 1:K
        y1(time_step, :) = y1(time_step, :) + coef(i) * exp(- eigen(i) * t(time_step)) * sin( i * pi * x );
    end
end

for time_step = 1:M
    for i = 1:K
        y2(time_step, :) = y2(time_step, :) + coef(i) * exp(- eigen(i) * t(time_step)) * cos( i * pi * x );
    end
end


[U1_, S1_, V1_] = svd(y1);
[U2_, S2_, V2_] = svd(y2);


%% 
for p = 1:8
    f_coef = zeros(1, K);
    f_coef(p) = 1;
    f_coef = f_coef / norm(f_coef);

    L = 1;
    tau = linspace(1e-1, 1e-1, L);

    z = zeros(L, N);
    z2 = zeros(L, N);
    eta = zeros(L, 1);
    eta2 = zeros(L, 1);

    for j = 1:L
        for i = 1:K
            z(j,:) = z(j,:) + f_coef(i) * exp(-eigen(i)* tau(j)) * cos(i * pi * x);
            z2(j,:) = z2(j,:) + f_coef(i) * exp(-eigen(i)* tau(j)) * sin(i * pi * x);
        end
    end

    for j=1:L
        eta(j) = ( norm( z(j,:) * V2_(:,1:26)) - norm(z(j,:)) ) / norm(z(j,:));
        eta2(j) = ( norm( z2(j,:) * V1_(:,1:26)) - norm(z2(j,:)) ) / norm(z2(j,:));
    end


    fprintf('p=%d, error=%e\n',...
        p, max(sqrt (abs(eta.^2 + abs(eta2.^2 ))  )));
end

