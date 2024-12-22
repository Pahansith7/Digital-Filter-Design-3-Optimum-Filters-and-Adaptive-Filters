% Adaptive Filter Implementation

clc;
clear all;

% Data Construction

% Generate Noise free signal (desired signal)
fs = 500;                             % Sampling frequency (Hz)
N = 5000;                             % Number of ponits
t = linspace(0, 10, N)';              % Time vector with fs =500 Hz

% Generate Noise free signal (desired signal)
yi = sawtooth(2*pi*2*t(1:N,1), 0.5);  % width = 0.5

% Generate non-stationary noise (n(n))
phi = pi/4;
n1 = 0.2*sin(2*pi*50*t(1:N/2,1)-phi);
n2 = 0.3*sin(2*pi*100*t(N/2+1:N,1)-phi);
n_non_stationary = [n1; n2];

% Gaussian white noise
snr = 20;
nwg = yi - awgn(yi, snr, 'measured');

n = n_non_stationary + nwg;         % Total noise
x = n + yi;                         % Primary input signal (n(n) + yi(n))

% Reference signal
a = 0.3;
phi1 = pi/6;
phi2 = pi/4;
r = a * (nwg + sin(2*pi*50*t + phi1) + sin(2*pi*100*t + phi2));

%--------------------------------------------------------------------------

% 2.1 (a)

% apply LMS filter

order_lms = 10;                   % order of LMS filter (arbitary value)
mu = 0.1;                         % step size (arbitary value)
[y_hat_lms, e_lms, w_lms] = lms_filter(x, r, order_lms, mu);

%--------------------------------------------------------------------------

% 2.1 (b)

figure(1);
subplot(4, 1, 1);
plot(t, yi);
title('Desired Signal y_i(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 2);
plot(t, x);
title('Input Signal x(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 3);
plot(t, e_lms);
title('LMS : Error Signal e(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 4);
plot(t, abs(yi - e_lms));
title('LMS : Absolute Error |e(n) - y(n)|');
xlabel('Samples');
ylabel('Absolute Error');

%--------------------------------------------------------------------------

% 2.1 (c)

M_values_lms = 1:1:20;             % Filter order range
mu_values = 0.001:0.001:0.1;       % Step size range (rate of convergence)

% Preallocate z matrix to store the MSE values
mse_values = zeros(length(mu_values), length(M_values_lms));


% Loop through both M_values_lms and mu_values to calculate MSE
for i = 1:length(M_values_lms)       
    for j = 1:length(mu_values)      
        order_lms = M_values_lms(i); 
        mu = mu_values(j);           
        
        % Call LMS filter function 
        [y_hat_lms, e_lms, w_lms] = lms_filter(x, r, order_lms, mu);
        
        % Compute the mean squared error (MSE) 
        mse_values(j, i) = mean((yi-e_lms).^2);  
    end
end


figure(2);
surf(M_values_lms, mu_values, mse_values);
xlabel('Filter Order (M)');
ylabel('Step size (mu)');
zlabel('Mean Squared Error (MSE)');
title('MSE vs Filter Order and Step Size');
grid on;
shading interp;  
lighting phong;

% Find the minimum value in the MSE matrix
minMSE = min(mse_values(:));

% Find the indices of the minimum value
[row, col] = find(mse_values == minMSE);

fprintf('Minimum MSE: %f at filter order = %f, mu = %f\n', minMSE, M_values_lms(col), mu_values(row));

% Apply optimum filter
opt_order_lms = M_values_lms(col);
opt_mu = mu_values(row);
[opt_y_hat_lms, opt_e_lms, opt_w_lms] = lms_filter(x, r, opt_order_lms, opt_mu);

figure(3);
subplot(4, 1, 1);
plot(t, yi);
title('Desired Signal y_i(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 2);
plot(t, x);
title('Input Signal x(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 3);
plot(t, opt_e_lms);
title('LMS : Error Signal e(n)for optimal filter');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 4);
plot(t, abs(yi - opt_e_lms));
title('LMS : Absolute Error |e(n) - y(n)|');
xlabel('Samples');
ylabel('Absolute Error');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply RLS filter

% 2.2 (a)

order_rls = 10;                         % order of RLS filter (arbitary value)
lambda = 0.9;                           % forgetting factor (arbitary value)
[y_hat_rls, e_rls, w_rls] = rls_filter(x, r, order_rls, lambda);

%--------------------------------------------------------------------------

% 2.2 (b)

figure(4);
subplot(4, 1, 1);
plot(t, yi);
title('Desired Signal y_i(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 2);
plot(t, x);
title('Input Signal x(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 3);
plot(t, e_rls);
title('RLS : Error Signal e(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

% Absolute error plot
subplot(4, 1, 4);
plot(t, abs(yi - e_rls));
title('RLS : Absolute Error |e(n) - y(n)|');
xlabel('Samples');
ylabel('Absolute Error');

%--------------------------------------------------------------------------

% 2.2 (c)

M_values_rls = 1:1:20;                      % Filter order range
lambda_values = 0.98:0.001:1;                 % forgetting factor

% Preallocate z matrix to store the MSE values
mse_values = zeros(length(lambda_values), length(M_values_rls));


for i = 1:length(M_values_rls)       
    for j = 1:length(lambda_values)      
        order_rls = M_values_rls(i); 
        lambda = lambda_values(j);          
        
        % Call LMS filter function 
        [y_hat_rls, e_rls, w_rls] = rls_filter(x, r, order_rls, lambda);
        
        % Compute the mean squared error (MSE) 
        mse_values(j, i) = mean((yi-e_rls).^2);  
    end
end


figure(5);
surf(M_values_rls, lambda_values, mse_values);
xlabel('Filter Order (M)');
ylabel('Step size (mu)');
zlabel('Mean Squared Error (MSE)');
title('MSE vs Filter Order and forgetting factor');
grid on;
shading interp;  
lighting phong;  


% Find the minimum value in the MSE matrix
minMSE = min(mse_values(:));

% Find the indices of the minimum value
[row, col] = find(mse_values == minMSE);

fprintf('Minimum MSE: %f at filter order = %f, lambda = %f\n', minMSE, M_values_rls(col), lambda_values(row));

% Apply optimum filter
opt_order_rls = M_values_rls(col);
opt_lambda = lambda_values(row); 
[opt_y_hat_rls, opt_e_rls, opt_w_rls] = rls_filter(x, r, opt_order_rls, opt_lambda);

figure(6);
subplot(4, 1, 1);
plot(t, yi);
title('Desired Signal y_i(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 2);
plot(t, x);
title('Input Signal x(n)');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

subplot(4, 1, 3);
plot(t, opt_e_rls);
title('RLS : Error Signal e(n) for optimu filter');
xlabel('Samples');
ylabel('Amplitude');
ylim([-2,2]);

% Absolute error plot
subplot(4, 1, 4);
plot(t, abs(yi - e_rls));
title('RLS : Absolute Error |e(n) - y(n)|');
xlabel('Samples');
ylabel('Absolute Error')

%--------------------------------------------------------------------------

% 2.2 (d)

load('idealECG')
yi  = idealECG(1:5000);                % Take 5000 sampleas as input signal
M = length(yi);                        % Length of the signal
Fs = 500;                              % sampling frequency(Hz)
t = (0:M-1)/Fs;                        % Time vector for the signal

x = yi' + n;                            % Noise remains same

% Apply LMS filter
[~, lms_filtered_ECG, ~] = lms_filter(x, r, opt_order_lms, opt_mu);

figure(7);
subplot(3, 1, 1);
plot(t(2250:2750), yi(2250:2750));
title('Ideal ECG yi(n)');
xlabel('Samples');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(t(2250:2750), x(2250:2750));
title('Input ECG Signal x(n)');
xlabel('Samples');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t(2250:2750) ,lms_filtered_ECG(2250:2750));
title('LMS : Filtered ECG');
xlabel('Samples');
ylabel('Amplitude');


% Apply RLS filter
[~, rls_filtered_ECG, ~] = rls_filter(x, r, opt_order_rls, opt_lambda);

figure(8);
subplot(3, 1, 1);
plot(t(2250:2750), yi(2250:2750));
title('Ideal ECG yi(n)');
xlabel('Samples');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(t(2250:2750), x(2250:2750));
title('Input ECG Signal x(n)');
xlabel('Samples');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t(2250:2750), rls_filtered_ECG(2250:2750));
title('RLS : Filtered ECG');
xlabel('Samples');
ylabel('Amplitude');



