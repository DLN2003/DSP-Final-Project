%% 4.2a) Same as 4.1a) for new filter where wc = 0.25*pi and L = 41. Also, measure response at w = {0, 0.1*pi, 0.25*pi, 0.4*pi, 0.5*pi, 0.75*pi}
clear

% Begin Filter analysis
L = 41; % Length of Filter
N = 2048; % # of Points for the DFT
w = -pi: 2*pi/N: pi - 2*pi/N; % Define frequency range 
wc = 0.25*pi; % Cutoff frequency for the filter

% Use Project function to define better bandpass filter 
H = BPFbetter(wc, L, N);
idx = N/2 + 1; % Define start index for positive frequencies

% Plot magnitude of the frequency response
figure; clf;
plot(w(idx:end), abs(H(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.25π ~ 0.7854
xline(0.7854, '--r', '\omega = 0.25\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% This Plot shows the positive index of our bandpass better filter with the main
% lobe centered at wc = 0.25*pi. We see very low side lobes coming off of our main 
% lobe. This suggests our passband threshold can be much closer to 0   

% Plot phase for frequency response 
figure; clf;
plot(w(idx:end), angle(H(idx:end)));
title('Phase of Frequency Response');
xlabel('\omega (rad)');
ylabel('\Theta(\omega)');

% This is the Phase of our passband which is linear in the region of each
% of the main lobe, but is not consistently linear in the regions of the
% side lobes. For our FIR filter this is as expected for a non normalized signal.



%% Use find function to calculate frequency response at the desired frequencies

% Define intermediate function for response at positive frequencies for
% indexing purposes
G = H(idx:end);

% Define values of omega
selected_w = [0, 0.1*pi, 0.25*pi, 0.4*pi, 0.5*pi, 0.75*pi];

% Initialize an array to store the magnitudes at selected frequencies
H_for_selected_w = zeros(1, length(selected_w));

% Initialize phase for part 4.2c)
Phase_for_selected_w = zeros(1, length(selected_w));

% Initalize an array of values for the target indicies to find H
target_idx = zeros(1, length(selected_w));

% Calculate H by finding the correct index in the range of frequencies
for i = 1:length(selected_w)
    % Find all the indices in w(idx:end) closest to the desired frequency
    target_idx = find(abs(w(idx:end) - selected_w(i)) <= 0.0013);
    
    % Getting the magnitude and phase values for the target indices
    H_for_selected_w(i) = abs(G(target_idx));
    Phase_for_selected_w(i) = angle(G(target_idx));
end

% Display the magnitudes for the selected frequencies
fprintf('|H(\x03c9)| for selected values of \x03c9:\n');
for i = 1:length(selected_w)
    fprintf('  %.2f\x03c0: %.4f\n', selected_w(i)/pi, H_for_selected_w(i));
end

% Display values of phase for part 4.2c)
fprintf('\x0398(\x03c9) for selected values of \x03c9:\n');
for i = 1:length(selected_w)
    fprintf('  %.2f\x03c0: %.4f\n', selected_w(i)/pi, Phase_for_selected_w(i));
end

% This list of Responses makes sense, you can just look at the graph. Most
% of the values of w are close to zero, but the value at wc = 0.25*pi is
% much higher (at 10.88) since that is the center of the passband and all the other
% selected frequencies fall outside the passband.
