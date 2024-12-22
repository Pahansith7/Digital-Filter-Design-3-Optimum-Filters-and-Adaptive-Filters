function [y_hat, e, w] = rls_filter(x, r, order, lambda)
    % RLS Adaptive Filter Function
    % Inputs:
    %   x      - Primary input signal (N x 1)
    %   r      - Reference signal (N x 1)
    %   M      - Filter order
    %   lambda - Forgetting factor (close to 1)
    % Outputs:
    %   y_hat  - Noise estimation (N x 1)
    %   e      - Error signal(filtered signal) (N x 1)
    %   w      - Final filter weights (M x 1)


    N = length(x);
    M = order + 1; % length of the filter

    % Initialize variables
    delta = 1e-2;             % Small constant for initialization of P(0)
    P = delta^(-1) * eye(M);  % Initial inverse correlation matrix P(0)
    w = zeros(M, 1);          % Initial weight vector w(0)
    e = zeros(N, 1);          % Error signal( filtered signal)
    y_hat = zeros(N, 1);      % Estimated noiose

    % Adaptive filtering process using RLS
    for n = M:N
        % Form the reference vector r_n (most recent M samples)
        r_n = r(n:-1:n-M+1);
        k = (lambda^(-1) * P * r_n) / (1 + lambda^(-1) * r_n' * P * r_n);
        P = lambda^(-1) * (P - k * r_n' * P);
        e(n) = x(n) - w' * r_n;
        w = w + k * e(n);
        y_hat(n) = w' * r_n;
    end
end
