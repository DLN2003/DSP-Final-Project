% Authored By: David Needens
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 5.3 Piano Octave Decoding
% 
% This MATLAB code sets up multiple bandpass filters based on predefined 
% frequency bands, calculates their frequency responses using a Hamming 
% window function, applies them to the signal xx in the time domain, 
% and stores the filtered outputs for each band.
%

%% 5.3b)
% Frequency bands for 5 BPF
bands = [BP_Filters(2, "StartingFreq_Hz_"), BP_Filters(2, "EndingFreq_Hz_"); 
        BP_Filters(3, "StartingFreq_Hz_"), BP_Filters(3, "EndingFreq_Hz_"); 
        BP_Filters(4, "StartingFreq_Hz_"), BP_Filters(4, "EndingFreq_Hz_");
        BP_Filters(5, "StartingFreq_Hz_"), BP_Filters(5, "EndingFreq_Hz_");
        BP_Filters(6, "StartingFreq_Hz_"), BP_Filters(6, "EndingFreq_Hz_");];
L = 101;                 % Filter length
N = 1024;                % FFT size for frequency response
numBands = size(bands, 1);
filterOutputs = zeros(numBands, length(xx)); % To store filter outputs

for i = 1:numBands
    % Calculate center frequency and normalized cutoff
    fc = (bands{i, "StartingFreq_Hz_"} + bands{i, "EndingFreq_Hz_"}) / 2;
    wc = 2 * pi * fc / fs; % Convert to normalized frequency
    
    % Call HammingNorm to generate frequency response
    H = HammingNorm(wc, L, N); % Obtain frequency response
    
    % Apply the filter (IFFT to create time-domain filter coefficients)
    h = ifft(ifftshift(H), 'symmetric');
    h = h(1:L); % Trim to filter length
    
    % Filter the signal
    filterOutputs(i, :) = filter(h, 1, xx);
end
