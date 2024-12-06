% Authored By: Giovanni Gutierrez
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises
%
% This code defines a function that creates a BPF 
% with a specified cutoff frequency, filter length, 
% and DFT length; using a Hamming Window to reduce
% noise from side lobes.
% 



%% Better Bandpass filter (Used first in 4.2a)

% Function for bandpass filter frequency response
function H = BPFbetter(wc, L, N)
    h = zeros(1, L); % Initialize impulse response
    
    % For loop defining impulse Response for BPF with Hamming window
    for n = 0:L-1
        h(n+1) = (0.54 - 0.46 * cos(2 * pi * n / (L - 1))) * cos(wc * (n - (L - 1) / 2)); 
    end

    % Frequency response for the bandpass filter
H = fftshift(fft(h,N)); % Take N-point DFT of Bandpass filter and shift 0 frequency to center

end
