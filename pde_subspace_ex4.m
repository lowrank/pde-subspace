N = 2^10;
M = 2^10;
K = 2000;

c= 2^0.125;

x = linspace(0, 1, N);
y = linspace(0, 1/c, N);
t = linspace(0, 1, M);

%% 
solx = zeros(M, N);
soly = zeros(M, N);
sol = zeros(N, N, M);
coefx = zeros(K, 1);
coefy = zeros(K, 1);
eigenx = zeros(K, 1);
eigeny = zeros(K, 1);

tic;
for i = 1:K
    coefx(i) = 1/(i^2);
    coefy(i) = 1/(i^2);
    eigenx(i) = pi^2 * (i^4);
    eigeny(i) = pi^2 * sqrt(2) * i^4;    
end
toc;

%%
tic;
for time_step = 1:M
    for i = 1:K
        solx(time_step, :) = ...
        solx(time_step, :) + ...
        coefx(i) * ...
        exp(- eigenx(i) * t(time_step)) * ...
        sin( i * pi * x );
    end
    
    for i = 1:K
        soly(time_step, :) = ...
            soly(time_step, :) + ...
            coefy(i) * ...
            exp(-eigeny(i) * t(time_step)) * ...
            sin(c * i * pi * y);
    end
end
toc;
%%
tic;
for time_step =1:M
    for x_id = 1:N
        for y_id=1:N
            sol(x_id, y_id, time_step) = solx(time_step, x_id) * soly(time_step, y_id);
        end
    end
end
toc;
%%
tic;
sol = reshape(sol, N*N, M);

[U_, S_, V_] = svds(sol, 600);
toc;

%% 
L = 1;
tau = linspace(1e-1, 1e-1, L);
small_z = zeros(L,1);

for p = 1:5
    for q = 1:5

        z = zeros(N, N, L);
        eta = zeros(L, 1);


        for j = 1:L
            for x_id = 1:N
                for y_id = 1:N
                    z(x_id, y_id, j) = ...
                        exp(-eigenx(p) * tau(j)) ...
                        * sin(p * pi * x(x_id)) * ...
                        exp(-eigeny(q) * tau(j)) * ...
                        sin(q * c * pi * y(y_id));
                end
            end
        end
        
        z = reshape(z, N*N, L);
        for j=1:L
            eta(j) = ( norm( z(:, j)' * U_(:,1:23)) - norm(z(:,j)) ) / norm(z(:, j));    
        end

        fprintf('p=%d, q=%d, error=%e\n', p, q, max(abs(eta)));
    end
end

