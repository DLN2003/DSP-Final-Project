%% Simple Band Pass Filter (used in 4.1a)

% Function for bandpass filter frequency response
function H = BPFsimp(wc, L, N)
    h = zeros(1, L); % Initialize impulse response
    % For loop defining impulse Response for BPF
    for n = 0:L-1
        h(n+1) = (2 / L) * cos(wc * n); 
    end

    % Frequency response for the bandpass filter
H = fftshift(fft(h,N)); % Take N-point DFT of Bandpass filter and shift 0 frequency to center

end
