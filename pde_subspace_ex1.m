N = 2^10;
M = 2^12;
K = 2000;

x = linspace(0, 1, N);
t = linspace(0, 1, M);
%% 
y = zeros(M, N);

coef = 1.0 ./ (1:K).^2 .* (-1).^(1:K);

for time_step = 1:M
    for i = 1:K
        y(time_step, :) = y(time_step, :) + coef(i) * exp(-pi^2 * i^2 * t(time_step)) * sin( i * pi * x );
    end
end

[U_, S_, V_] = svd(y);

for p = 1:8
    f_coef = zeros(1, K);
    f_coef(p) = 1;
    f_coef = f_coef / norm(f_coef);

    L = 1;
    tau = linspace(1e-1, 1e-1, L);

    z = zeros(L, N);
    eta = zeros(L, 1);

    for j = 1:L
        for i = 1:K
            z(j,:) = z(j,:) + f_coef(i) * exp(-pi^2 * i^2 * tau(j)) * sin(i * pi * x);
        end
    end

    for j=1:L
        eta(j) = ( norm( z(j,:) * V_(:,1:20)) - norm(z(j,:)) ) / norm(z(j,:));
    end

    fprintf('p=%d, error=%e\n',...
        p, max(abs(eta)));

end

