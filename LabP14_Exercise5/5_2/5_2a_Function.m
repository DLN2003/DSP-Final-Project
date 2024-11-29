%% 5.2a) 
% Normalized Hanning Bandpass Filter (5.2)

% Same as function "BPFbetter" but with a scaling term "B" so that the
% maximum value at the center frequency is equal to one.
function H = HammingNorm(wc, L, N)
    h = zeros(1, L); % Initialize impulse response
    
    % For loop defining impulse Response for BPF with Hamming window
    for n = 0:L-1
        h(n+1) = (0.54 - 0.46 * cos(2 * pi * n / (L - 1))) * cos(wc * (n - (L - 1) / 2)); 
    end
    % Frequency response for the bandpass filter
H1 = fftshift(fft(h,N)); % Take N-point DFT of Bandpass filter and shift 0 frequency to center
B = 1/max(abs(H1));% Scaling factor so that max(abs(H)) = 1.
H = B*H1;
