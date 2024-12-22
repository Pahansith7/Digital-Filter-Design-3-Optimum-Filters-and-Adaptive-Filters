function w = computeWienerWeights(yi, n, M)
    % Compute wiener filter coefficients 
    % Inputs:
    %   y_i - Input signal (N x 1)
    %   n   - Noise signal (n x 1)
    %   M   - Filter order
    % Outputs:
    %   w   - Filter weights (M x 1)
 
    % Compute the auto-correlation of the input signal
    R_YY = xcorr(yi, M-1, 'unbiased');
    R_YY_matrix = toeplitz(R_YY(M:end)); 

    % Compute the auto-correlation of the noise signal
    R_NN = xcorr(n, M-1, 'unbiased');   
    R_NN_matrix = toeplitz(R_NN(M:end)); 

    % Compute the cross-correlation
    R_Yy = R_YY(M:end); 

    % Compute the Wiener filter weights
    w = inv(R_YY_matrix + R_NN_matrix) * R_Yy.';
    

end
