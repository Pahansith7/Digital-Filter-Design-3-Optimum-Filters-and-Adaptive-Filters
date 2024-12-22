function [y_hat, e, w] = lms_filter(x, r, order, mu)
    % LMS Adaptive Filter Function
    % Inputs:
    %   x  - Primary input signal (N x 1)
    %   r  - Reference signal (N x 1)
    %   M  - Filter order
    %   mu - Step size (learning rate)
    % Outputs:
    %   y_hat - Noise estimation (N x 1)
    %   e     - Error signal (N x 1)
    %   w     - Filter weights (M x 1)

    N = length(x);
    M = order+1;

    % Initialize variables
    w = zeros(M, 1);  % Filter weights initialized to zero
    e = zeros(N, 1);  % Estimated signal
    y_hat = zeros(N, 1);  % Estimated noise

    for n = M:N
        r_n = r(n:-1:n-M+1);
        y_hat(n) = w' * r_n;
        e(n) = x(n) - y_hat(n);
        w = w + 2 * mu * e(n) * r_n;
    end
end
