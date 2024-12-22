function [yhat_f, W_f] = wiener_filter_fd(x, yi, n)
    %  Compute wiener filter coefficients and Wiener filterd signal in
    %  frequency domain
    % Inputs:
    %   x   - Primary input signal
    %   yi  - Input signal (N x 1)
    %   n   - Noise signal (n x 1)
    %   M   - Filter order
    % Outputs:
    %   w_f   - Filter weights (M x 1)
    %   yhat_f   - Filtered signal (M x 1)
    
    % Get the length of the input signal 
    len_x = length(x);

    S_yy = abs(fft(yi, 2*len_x-1)).^2;  
    S_nn = abs(fft(n, 2*len_x-1)).^2;   

    % Wiener filter calculation in frequency domain
    W_f = S_yy ./ (S_yy + S_nn); 

    % Apply the Wiener filter to the input signal in frequency domain
    S_xx = fft(x, 2*len_x-1);  
    S_hat_yy = W_f .* S_xx;    

    % Inverse FFT to get the filtered signal back in the time domain
    yhat_f = ifft(S_hat_yy);
    yhat_f = yhat_f(1:len_x); 

end
