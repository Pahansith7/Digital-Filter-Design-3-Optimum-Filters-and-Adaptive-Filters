% 1. Wiener filter

clear all;
clc;

% Data Construction
load('idealECG.mat');
yi = idealECG;                              % ideal ECG signal
M = length(yi);                             % Length of the signal
Fs = 500;                                   % sampling frequency(Hz)
t = (0:M-1)/Fs;                             % Time vector for the signal

SNR = 10;
rng(0);
n_wg = awgn(yi, SNR, 'measured') - yi;      % white Gaussian noise

n_50 = 0.2 * sin(2 * pi * 50 * t);          % 50Hz noise

n  = n_wg + n_50;                           % Total noise
x = yi + n;                                 % Primary input signal

% Plot the constructed signals
figure(1);
subplot(4,1,1);
plot(t(1:500), yi(1:500));
title('Ideal ECG Signal yi(n)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,2);
plot(t(1:500), n_wg(1:500));
title('White Gaussian Noise eta_{wg}(n)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,3);
plot(t(1:500), n_50(1:500));
title('50 Hz Sinusoidal Noise eta_{50}(n)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,4);
plot(t(1:500), x(1:500));
title('Noisy Input Signal x(n) = yi(n) + eta_{wg}(n) + eta_{50}(n)');
xlabel('Time (s)');
ylabel('Amplitude');

%--------------------------------------------------------------------------

% Discrete time-domain implementation of the Wiener filter

% Part 1

% 1.1 (a)

M = 20;                                % Filter Length (an arbitary value)

yi_1 = yi(26:145);                     % Desired signal
yi_1 = yi_1 - mean(yi_1);              % For E[y] to be zero.
figure(2)
plot(yi_1);
title('Desired Signal y(n) - A single ECG Beat ');
xlabel('Time (s)');
ylabel('Amplitude');

n_ = x(27:46);                         % Isoelectric segment
n_1 =  [n_ n_ n_ n_ n_ n_];            % Noise signal
n_1 = (n_1- mean(n_1));
figure(3)
plot(n_1);
title('Noise Signal n(n)');
xlabel('Time (s)');
ylabel('Amplitude');

w_1 = computeWienerWeights(yi_1, n_1, M);

x_1 = x - mean(x);                    % Input Signal
x_1 = x_1 - mean(x_1);

yhat_1 = filter(w_1, 1, x_1);

figure(4)
plot(t(1:500), x_1(1:500), 'g', 'LineWidth', 2);    
hold on; 
plot(t(1:500), yhat_1(1:500), 'r', 'LineWidth', 2);  
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2);  
hold off;  
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Wiener Filtered signal for an arbitary order')


% % Find auto-correlation using x(n) and yi(n)

% R_XX = xcorr(x, M-1);   
% R_XX_matrix = toeplitz(R_XX(M:end));
% 
% % Compute the cross-correlation between input and desired signal
% R_Xy = xcorr(yi, x, M-1);  
% R_Xy = R_Xy(M:end);    
% 
% % Solve for the Wiener filter coefficients w
% w = inv(R_XX_matrix)*R_Xy.'  % Solve Wiener-Hopf equation
% plot(filter(w, 1, x))  
% plot(x) 

%--------------------------------------------------------------------------

% 1.1 (b)

% Obtain the optimum filter order

errors = [];  % Initialize an empty array to store RMS errors

for M_i = 1:120
    
    w_i = computeWienerWeights(yi_1, n_1, M_i);
    
    signal = filter(w_i, 1, x_1); 

    error = yi - signal;  
    
    rms_error = sqrt(mean(error.^2));  

    errors = [errors, rms_error];  
end

figure(5);
plot(1:120, errors);
title('RMS Error vs M_i');
xlabel('M_i');
ylabel('RMS Error');
grid on;

M_optimum = 41;                              % Obtained Using the grapgh
w_optimum = computeWienerWeights(yi_1, n_1, M_optimum);

% fvtool(w_optimum, 1)

%--------------------------------------------------------------------------

% 1.1 (c)

% FIlter the noisy signal using the optimum filter

yhat_optimum_1 = filter(w_optimum, 1, x_1); % Apply the filter

figure(6)
plot(t(1:500), x_1(1:500), 'g', 'LineWidth', 2);  
hold on;  
plot(t(1:500), yhat_optimum_1(1:500), 'r', 'LineWidth', 2);
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2); 
hold off;  
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Wiener Filtered signal for optimum filter order')

%--------------------------------------------------------------------------

% 1.1 (d)

% Plot each signal's PSDs

figure(7);
[p_yi, f] = periodogram(yi, [], [], Fs);
[p_n, ~] = periodogram(n, [], [], Fs);
[p_x, ~] = periodogram(x, [], [], Fs);
[p_yhat_optimum, ~] = periodogram(yhat_optimum_1, [], [], Fs);
plot(f, 10*log10(p_yi), f, 10*log10(p_n), f, 10*log10(p_x), f, 10*log10(p_yhat_optimum)); 
legend('yi(n)', 'n(n)', 'x(n)', 'yhat(n)');
title('Power Spectral Density of Signals');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% part 2

% Construct a mathematical ECG model
len = 120;                           % Total length of the ECG signal
R_peak_height = 1.5;                 % Height of the R peak
P_wave_height = 0.1;                 % Height of the P wave
T_wave_height = 0.25;                % Height of the T wave

liner_ecg_signal = zeros(1, len);   % Initialize the ECG signal

% Generate P Wave
p_wave = P_wave_height * (1 - cos(linspace(0, 2*pi, 15))); 
liner_ecg_signal(20:34) = liner_ecg_signal(20:34) + p_wave;

% Generate QRS Complex 
q_wave = linspace(0, -0.2, 5);         % Q wave (small downward slope)
r_wave = [linspace(-0.2, R_peak_height, 2), linspace(R_peak_height, -0.3, 3)]; % Triangular R wave
s_wave = linspace(-0.3, 0, 10);        % S wave (downward slope)

% Combine QRS components into the signal
liner_ecg_signal(40:44) = liner_ecg_signal(40:44) + q_wave;  % Add Q wave
liner_ecg_signal(44:48) = liner_ecg_signal(44:48) + r_wave;  % Add R wave 
liner_ecg_signal(48:57) = liner_ecg_signal(48:57) + s_wave;  % Add S wave

% Generate T Wave 
t_wave = T_wave_height * (1 - cos(linspace(0, 2*pi, 21)));
liner_ecg_signal(80:100) = liner_ecg_signal(80:100) + t_wave;

% 1.1 (a)

M = 20;                                 % Filter Length (an arbitary value)

yi_2 = liner_ecg_signal;                % Desired signal

yi_2 = yi_2 - mean(yi_2);               % For E[y] to be zero

figure(8)
plot(yi_2);
title('Desired Signal y(n) - A single ECG Beat ');
xlabel('Time (s)');
ylabel('Amplitude');

n_ = x(27:46);                          % Isoelectric segment
n_2 =  [n_ n_ n_ n_ n_ n_];             % Noise signal
n_2 = (n_2- mean(n_2));
figure(9)
plot(n_2);
title('Noise Signal n(n)');
xlabel('Time (s)');
ylabel('Amplitude');

w_2 = computeWienerWeights(yi_2, n_2, M);

x_2 = x;
x_2 = x_2 - mean(x_2);

yhat_2 = filter(w_2, 1, x_2);

figure(10)
plot(t(1:500), x_2(1:500), 'g', 'LineWidth', 2);  
hold on;  
plot(t(1:500), yhat_2(1:500), 'r', 'LineWidth', 2);  
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2); 
hold off;
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Wiener Filtered signal for an arbitary order')

%--------------------------------------------------------------------------

% 1.1 (b)

errors = [];  

for M_i = 1:120
    w_i = computeWienerWeights(yi_2, n_2, M_i);
    
    signal = filter(w_i, 1, x_2); 
    
    error = yi - signal;  

    rms_error = sqrt(mean(error.^2));  
    
    errors = [errors, rms_error];  
end

figure(11);
plot(1:120, errors);
title('RMS Error vs M_i');
xlabel('M_i');
ylabel('RMS Error');
grid on;

M_optimum_2 = 76;
w_optimum_2 = computeWienerWeights(yi_2, n_1, M_optimum_2);
% fvtool(w_optimum_2, 1)

%--------------------------------------------------------------------------

% 1.1 (c)

yhat_optimum_2 = filter(w_optimum_2, 1, x_2);

figure(12)
plot(t(1:500), x_2(1:500), 'g', 'LineWidth', 2); 
hold on; 
plot(t(1:500), yhat_optimum_2(1:500), 'r', 'LineWidth', 2);
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2); 
hold off; 
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Wiener Filtered signal for optimum filter order')

%--------------------------------------------------------------------------

% 1.1 (d)

% Plot each signal's PSDs
figure(13);
[p_yi, f] = periodogram(yi, [], [], Fs);
[p_n, ~] = periodogram(n, [], [], Fs);
[p_x, ~] = periodogram(x, [], [], Fs);
[p_yhat_optimum_2, ~] = periodogram(yhat_optimum_2, [], [], Fs);
plot(f, 10*log10(p_yi), f, 10*log10(p_n), f, 10*log10(p_x), f, 10*log10(p_yhat_optimum_2)); 
legend('yi(n)', 'n(n)', 'x(n)', 'yhat(n)');
title('Power Spectral Density of Signals');
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Frequency domain implementation of the Wiener filter

% 1.2 (a)

[yhat_f_1, W_f_1] = wiener_filter_fd(x_1, yi_1, n_1);

figure(14)
plot(t(1:500), x_1(1:500), 'g', 'LineWidth', 2);
hold on; 
plot(t(1:500), yhat_f_1(1:500), 'r', 'LineWidth', 2); 
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2);
hold off; 
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Signal using Wiener Filter Derived in frequency Domain');

[yhat_f_2, w_f_2] = wiener_filter_fd(x_2, yi_2, n_2);

figure(15)
plot(t(1:500), x_2(1:500), 'g', 'LineWidth', 2);
hold on; 
plot(t(1:500), yhat_f_2(1:500), 'r', 'LineWidth', 2);
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2); 
hold off;
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Signal using Wiener Filter Derived in frequency Domain');

%--------------------------------------------------------------------------

% 1.2 (b)

% Calculate the Mean Squared Error (MSE)
MSE_1 = mean((yi - yhat_f_1).^2);
MSE_2 = mean((yi - yhat_f_2).^2);

figure(16)
plot(t(1:500), x_2(1:500), 'g', 'LineWidth', 2);
hold on; 
plot(t(1:500), yhat_f_1(1:500), 'r', 'LineWidth', 2);
plot(t(1:500), yhat_f_2(1:500), 'm', 'LineWidth', 2);
plot(t(1:500), yi(1:500), 'b', 'LineWidth', 2)
hold off;
legend(['Input signal x(i)'], ...
       ['Filtered signal - Part 1: MSE ' num2str(MSE_1)], ...
       ['Filtered signal - Part 2: MSE ' num2str(MSE_2)], ...
       ['Desired signal y(i)'], 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Signal using Wiener Filter Derved in frequency Domain');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Effect on non-stationary noise on Wiener filtering

N = length(yi);
n1 = 0.2*sin(2*pi*50*t(1:N/2));
n2 = 0.3*sin(2*pi*100*t(N/2+1:N));
n_non_st = [n1'; n2'];
x_non_st = yi + n_non_st';
x_non_st = x_non_st - mean(x_non_st);

%--------------------------------------------------------------------------

% 1.3 (a)

yhat_non_st_1 = filter(w_1, 1, x_non_st);

len_x_non_st = length(x_non_st);
S_xx_non_st = fft(x_non_st, 2*len_x_non_st-1); 
Shat_yy_non_st = W_f_1 .* S_xx_non_st;    
yhat_non_st_1_f = ifft(Shat_yy_non_st);
yhat_non_st_1_f = yhat_non_st_1_f(1:len_x_non_st);  

%--------------------------------------------------------------------------

% 1.3 (b)

figure(17)
plot(t(3590:4090), x_non_st(3590:4090), 'g', 'LineWidth', 2);
hold on; 
plot(t(3590:4090), yhat_non_st_1(3590:4090), 'r', 'LineWidth', 2);
plot(t(3590:4090), yi(3590:4090), 'b', 'LineWidth', 2); 
hold off;
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');

figure(18)
plot(t(3590:4090), x_non_st(3590:4090), 'g', 'LineWidth', 2);
hold on; 
plot(t(3590:4090), yhat_non_st_1_f(3590:4090), 'r', 'LineWidth', 2);
plot(t(3590:4090), yi(3590:4090), 'b', 'LineWidth', 2); 
hold off;
legend('Input signal x(i)', 'Filtered signal', 'Desired signal y(i)', 'Location', 'best');
xlabel('Time (s)');
ylabel('Amplitude');


