%% 5.2a) See function "HanningNorm(wc, L, N);"
%% 5.2b)
wc = BP_Filters.CenterFreq_omegaC_;
N = 2048; % # of Points for the DFT
w = -pi: 2*pi/N: pi - 2*pi/N;%Frequency range
% Filters for octaves 1-7. L was determined through trial and error and was
% chosen based on which L generated a bandwidth closest to the design
% bandwidth for each octave.
% L = [519, 251, 128,63,32,17,9];% Values originally found for 0.5 stopband
L = [943, 546, 246,133, 67, 34, 17]; % Values with a 0.01 stopband


%% 5.2c)
HOct7 = HammingNorm(wc(7),L(7),N);
HOct6 = HammingNorm(wc(6),L(6),N);
HOct5 = HammingNorm(wc(5),L(5),N);
HOct4 = HammingNorm(wc(4),L(4),N);
HOct3 = HammingNorm(wc(3),L(3),N);
HOct2 = HammingNorm(wc(2),L(2),N);
HOct1 = HammingNorm(wc(1),L(1),N);
HOctTot = [HOct1;HOct2;HOct3;HOct4;HOct5;HOct6;HOct7];
whalf = w(length(w)/2:end);
HOctHalf = HOctTot(:,length(w)/2:end);
figure(9)
plot(whalf,abs(HOctHalf),LineWidth=2)
hold on
xline(wc, '--r', strcat('\omega = ' , string(wc)), 'LabelOrientation', 'aligned', 'LabelVerticalAlignment','middle');

%% Function

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
end
