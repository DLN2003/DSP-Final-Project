%% Octave Band Filtering: Lab P-14: 4 Lab Exercises
%% 4.1 Simple Bandpass Filter 
% Use the impulse response of FIR: h(n) = (2/L)*cos(wc*n), 0 <= n < L
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

%% 4.1b) Using the 0.5 level passband, find width of passband 

% Use Project function to define the pass band width
Passband = PBWidth(H(idx:end),w(idx:end),0.5);

% Display passband width
fprintf('Passband width @ the 0.5 level is approximately %.4f radians for the L = 40 filter\n', Passband);

% We can imagine this width drawn on our first plot, in 4.1a), between the sides of
% the main lobe that is centered at wc. It can give us a range of
% frequencies accepted by our passband using the cutoff frequency where
% Range of accepted w = wc +- PassBand_Width/2 

%% 4.1c) Make plots to measure the passband at L = 20 and L = 80 with wc the same; measure their passband widths

% BPF for L = 20 using project function 
H20 = BPFsimp(wc, 20, N);

% Plot magnitude of the frequency response for L = 20 BPF
figure; clf;
plot(w(idx:end), abs(H20(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.4π ~ 1.257 
xline(1.257, '--r', '\omega = 0.4\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
yline(1, '--')

% Use Project function to define the pass band width
Passband = PBWidth(H20(idx:end),w(idx:end),0.5);

% Display passband width for L = 20 BPF
fprintf('Passband width @ the 0.5 level is approximately %.4f radians for the L = 20 filter\n', Passband);



%% BPF for L = 80 using project function 
H80 = BPFsimp(wc, 80, N);

% Plot magnitude of the frequency response for L = 80 BPF
figure; clf;
plot(w(idx:end), abs(H80(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.4π ~ 1.257 
xline(1.257, '--r', '\omega = 0.4\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
yline(1, '--')

% Use Project function to define the pass band width
Passband = PBWidth(H80(idx:end),w(idx:end),0.5);

% Display passband width for L = 80 BPF
fprintf('Passband width @ the 0.5 level is approximately %.4f radians for the L = 80 filter\n', Passband);

% It is obvious from the figures that as L increases, the width of the pass
% band narrows. At half the original length, the passband is twice as wide;
% at twice the original length, the passband is half as wide.

%% 4.2) A Better BPF
% Use a Hamming window to adjust the BPF where now:
%  h(n) = (0.54-0.46*cos(2*pi*n/(L-1)))cos(wc(n-(L-1)/2)), n = 1,2,3,...,L-1

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

%% 4.2b) Find the passband width at a threshold level of 50%. Then plot the graphs for two more BPFs at L = 21 and L = 81 and get their passbands.

% 50% of the maximum 10.88 is 5.44, so this defines our threshold. 

% This threshold can be calculated by
th = max(abs(G)) * 0.5;

% We can use the same function defined in part 4.1(b) 
Passband = PBWidth(G,w(idx:end),th);

% Display passband width for L = 41
fprintf('Passband width @ the 50%% level is approximately %.4f radians for the L = 41 filter\n', Passband);



%% We now define the BPF for L = 21
H21 = BPFbetter(wc, 21, N);
idx = N/2 + 1; % Define start index for positive frequencies

% Plot magnitude of the frequency response for L = 21
figure; clf;
plot(w(idx:end), abs(H21(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.25π ~ 0.7854
xline(0.7854, '--r', '\omega = 0.25\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% Set a threshold based on the peak value for H21
th21 = max(abs(H21(idx:end))) * 0.5;

% Use PBWidth with the updated threshold
Passband21 = PBWidth(H21(idx:end), w(idx:end), th21);

% Display passband width for L = 21
fprintf('Passband width @ the 50%% level is approximately %.4f radians for the L = 21 filter\n', Passband21);



%% And now for the analysis of the L = 81 BPF
H81 = BPFbetter(wc, 81, N);
idx = N/2 + 1; % Define start index for positive frequencies

% Plot magnitude of the frequency response for L = 81
figure; clf;
plot(w(idx:end), abs(H81(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.25π ~ 0.7854
xline(0.7854, '--r', '\omega = 0.25\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% Set a threshold based on the peak value for H81
th81 = max(abs(H81(idx:end))) * 0.5;

% Use PBWidth with the updated threshold
Passband81 = PBWidth(H81(idx:end), w(idx:end), th81);

% Display passband width for L = 81
fprintf('Passband width @ the 50%% level is approximately %.4f radians for the L = 81 filter\n', Passband81);

% The relationship is similar to that found in part 7.1 (c) in terms of the
% passband width. However, this filter is not normalized, so the magnitude
% of the frequency response changes depending on the value of L that you
% chose. Lower values of L will mean that the passband has less gain, and
% a higher value of L will correlate to more gain.


%% 4.2c) Given a specific input, determine the output signal by hand

% Display image of hand written derivation of output as requested by the problem.
filename = "Problem_4_2c.png";
imshow(filename)
%% 4.2d) Use frequency response to explain why the filter only passes at cutoff

% Observing the frequency response of the 41-length filter from 4.2a) we
% can see that the response is only significant at the 0.25*pi frequency.
% The response dies down to very shallow side lobes outside of this
% passband. If you observe the output equation from 4.2c) it is easy to see
% that the magnitude of the response for an input centered anywhere in the
% frequency range of the filter will be either amplified of attentuated by
% the filter. If an input component falls within the passband, it will
% be dominant in the output. However, if it falls within the stop-band, it
% will be greatly attenuated and not be greatly represented in the output.


%%  Project Functions

%% Simple Band Pass Filter (4.1a)

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


%% Width of pass band (4.1b)

% Function that finds passband width using the magnitude response, omega,
% and the threshold level for the passband.

function Passband = PBWidth(H,w,th)

% Find frequencies where |H(w)| is above threshold
Hmag = abs(H); % Defintion of magnitude response
H1 = find(Hmag >= th, 1, 'first'); % Find first index where |H(w)| is close to threshold
H2 = find(Hmag >= th, 1, 'last'); % Find last index where |H(w)| is close to threshold

% Convert indices to frequency values
w1 = w(H1); % Frequency at start of passband
w2 = w(H2); % Frequency at end of passband
Passband = w2 - w1; % Width of the passband
end

%% Better Bandpass filter (4.2a)

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