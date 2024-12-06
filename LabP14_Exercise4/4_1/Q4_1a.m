
% Authored By: Giovanni Gutierrez
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 4.1 Bandpass Filter Design
%
% To Decode the piano Octaves we need to create a filter.
% Bandpass filtering (BPF) is going to be the most useful
% in this process. This code exemplifies the use of a
% simple BPF to help understand the process involved
%


%% 4.1a) Generate Bandpass with wc = 0.4*pi and L = 40, plot magnitude and phase response


% Begin Filter analysis
clear
L = 40; % Length of Filter
N = 2048; % # of Points for the DFT
w = -pi: 2*pi/N: pi - 2*pi/N; % Define frequency range 
wc = 0.4*pi; % Cutoff frequency for the filter

% Frequency response for the bandpass filter
H = BPFsimp(wc, L, N); % Obtain frequency response of Bandpass filter
idx = N/2 + 1; % Define start index for positive frequencies

% Plot magnitude of the frequency response
figure; clf;
plot(w(idx:end), abs(H(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.4π ~ 1.257
xline(1.257, '--r', '\omega = 0.4\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
yline(1, '--')

% This Plot shows the positive index of our bandpass filter with the main
% lobe centered at wc. It's side lobes are reasonably large, meaning that
% the threshold level we choose for our passband cannot be too close to 0.

% Plot phase for frequency response 
figure; clf;
plot(w(idx:end), angle(H(idx:end)));
title('Phase of Frequency Response');
xlabel('\omega (rad)');
ylabel('\Theta(\omega)');

% This is the Phase of our passband which is linear in the region of each
% lobe, though the line centered at the cutoff frequency is larger &
% actually representative of our phase. This is expected for an FIR filter.
