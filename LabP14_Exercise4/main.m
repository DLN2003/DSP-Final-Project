%% Octave Band Filtering: Lab P-14: 4 Lab Exercises
clear;
close all;
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
figure(1); clf;
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
figure(2); clf;
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
figure(3); clf;
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



% BPF for L = 80 using project function 
H80 = BPFsimp(wc, 80, N);

% Plot magnitude of the frequency response for L = 80 BPF
figure(4); clf;
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
figure(5); clf;
plot(w(idx:end), abs(H(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.25π ~ 0.7854
xline(0.7854, '--r', '\omega = 0.4\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% This Plot shows the positive index of our bandpass better filter with the main
% lobe centered at wc = 0.25*pi. We see very low side lobes coming off of our main 
% lobe. This suggests our passband threshold can be much closer to 0   

% Plot phase for frequency response 
figure(6); clf;
plot(w(idx:end), angle(H(idx:end)));
title('Phase of Frequency Response');
xlabel('\omega (rad)');
ylabel('\Theta(\omega)');

% This is the Phase of our passband which is linear in the region of each
% of the main lobe, but is not consistently linear in the regions of the
% side lobes. For our FIR filter this is as expected for a non normalized signal.



% Use find function to calculate frequency response at the desired frequencies

% Define intermediate function for response at positive frequencies for
% indexing purposes
G = H(idx:end);

% Define values of omega
selected_w = [0, 0.1*pi, 0.25*pi, 0.4*pi, 0.5*pi, 0.75*pi];

% Initialize an array to store the magnitudes at selected frequencies
H_for_selected_w = zeros(1, length(selected_w));


% Initalize an array of values for the target indicies to find H
target_idx = zeros(1, length(selected_w));

% Calculate H by finding the correct index in the range of frequencies
for i = 1:length(selected_w)
    % Find all the indices in w(idx:end) closest to the desired frequency
    target_idx = find(abs(w(idx:end) - selected_w(i)) <= 1e-2);
    
    % Getting the H vlaues using our intermediate function G evaluated at
    % the fourth of the 6 indices found
    H_for_selected_w(i) = abs(G(target_idx(4)));
end

% Display the magnitudes for the selected frequencies
fprintf('|H(\x03c9)| for selected values of \x03c9:\n');
for i = 1:length(selected_w)
    fprintf('  %.2f\x03c0: %.4f\n', selected_w(i)/pi, H_for_selected_w(i));
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



% We now define the BPF for L = 21
H21 = BPFbetter(wc, 21, N);
idx = N/2 + 1; % Define start index for positive frequencies

% Plot magnitude of the frequency response for L = 21
figure(7); clf;
plot(w(idx:end), abs(H21(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.25π ~ 0.7854
xline(0.7854, '--r', '\omega = 0.4\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

% Set a threshold based on the peak value for H21
th21 = max(abs(H21(idx:end))) * 0.5;

% Use PBWidth with the updated threshold
Passband21 = PBWidth(H21(idx:end), w(idx:end), th21);

% Display passband width for L = 21
fprintf('Passband width @ the 50%% level is approximately %.4f radians for the L = 21 filter\n', Passband21);



% And now for the analysis of the L = 81 BPF
H81 = BPFbetter(wc, 81, N);
idx = N/2 + 1; % Define start index for positive frequencies

% Plot magnitude of the frequency response for L = 81
figure(8); clf;
plot(w(idx:end), abs(H81(idx:end)));
title('Magnitude of Frequency Response ');
xlabel('\omega (rad)');
ylabel('|H(\omega)|');

% Add a vertical line at ω = 0.25π ~ 0.7854
xline(0.7854, '--r', '\omega = 0.4\pi', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

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


%% 4.2c)

%% 4.2d)

%% 5.1
BP_Filters = readtable("Bandpass_Filters.xlsx"); % Load Bandpass filters from file
BP_Filters(8:11,:) = []; % Removes unnecessary rows
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

%% 5.3a)
fs = 8000;              % Sampling frequency
t = 0:1/fs:0.85;        % Time vector from 0 to 0.85 seconds with 1/fs step
xx = zeros(size(t));    % Initialize signal vector

% Define the three time intervals
interval1 = (t >= 0) & (t < 0.25);
interval2 = (t >= 0.3) & (t < 0.55);
interval3 = (t >= 0.6) & (t < 0.85);

% Add sinusoids in each time interval
xx(interval1) = cos(2*pi*220*t(interval1));
xx(interval2) = cos(2*pi*880*t(interval2));
xx(interval3) = cos(2*pi*440*t(interval3)) + cos(2*pi*1760*t(interval3));

% Plot the signal
figure(10);
plot(t, xx);
title('Generated Signal x(t)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

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
    
%% 5.3c)
% Plot Frequency Responses for All Filters
figure(11);
for i = 1:numBands
    subplot(numBands, 1, i);
    plot(t, filterOutputs(i, :));
    title(['Filter Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
        '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

%% 5.3d) Validate Output Signals by Comparing Magnitudes and Phases
% Frequency responses for validation
filterOutputMaxFreq = zeros(3, numBands);
regionIndex = [1, 0.25*fs; 0.3*fs, 0.55*fs; 0.6*fs, 0.85*fs];
regions = {xx(1:0.25*fs), xx(0.3*fs:0.55*fs), xx(0.6*fs:0.85*fs)}; % Cell array to store regions
regionOctave = zeros(3, 1); % To store octave results for each region

%For the input signal xx, based on the given definition:

%xx1 is a sinusoid at 220 Hz → Expected Octave = 3.
%xx2 is a sinusoid at 880 Hz → Expected Octave = 5.
%xx3 contains 440 Hz and 1760 Hz → Expected Octaves = 4.

% Validation of Regions
expectedOctaves = [3, 5, 4]; % Expected octaves for xx1, xx2, xx3

% Define frequency axis limits
f = (0:N-1) * (fs / N);
maxFreqIndex = find(f <= BP_Filters{6, "EndingFreq_Hz_"}); 
maxFreqIndex = maxFreqIndex(end);
minFreqIndex = find(f >= BP_Filters{1, "StartingFreq_Hz_"}); 
minFreqIndex = minFreqIndex(1);

for regionIdx = 1:3
    currentRegion = regions{regionIdx}; % Select the current region (xx1, xx2, xx3)
    
    % Initialize storage for filter outputs in the current region
    filterOutputs = zeros(numBands, length(currentRegion));
    maxMagnitudes = zeros(1, numBands); % To store the max magnitude for each filter band
    
    % Filter each region
    for i = 1:numBands
        % Calculate center frequency and normalized cutoff
        fc = (bands{i, "StartingFreq_Hz_"} + bands{i, "EndingFreq_Hz_"}) / 2;
        wc = 2 * pi * fc / fs; % Convert to normalized frequency

        % Call HammingNorm to generate frequency response
        H = HammingNorm(wc, L, N); % Obtain frequency response

        % Apply the filter (IFFT to create time-domain filter coefficients)
        h = ifft(ifftshift(H), 'symmetric');
        h = h(1:L); % Trim to filter length

        % Filter the current region signal
        filterOutputs(i, :) = filter(h, 1, currentRegion);
    end

    % Compute and store the max frequencies for validation
    figure(regionIdx + 11); % Start figure numbering from 12
    sgtitle(['Region ', num2str(regionIdx)]);

    for i = 1:numBands
        % Fourier Transform of the filter output
        Y = fft(filterOutputs(i, :), N);
        Y = Y(minFreqIndex:maxFreqIndex);
        f_band = f(minFreqIndex:maxFreqIndex);

        % Magnitude and phase plots
        subplot(numBands, 2, 2*i-1);
        plot(f_band, abs(Y));
        title(['Magnitude of Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
            '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
        xlabel('Frequency (Hz)');
        xlim([f_band(1) f_band(end)]);
        ylabel('Magnitude');
        grid on;

        subplot(numBands, 2, 2*i);
        plot(f_band, angle(Y) / pi);
        title(['Phase of Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
            '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
        xlabel('Frequency (Hz)');
        xlim([f_band(1) f_band(end)]);
        ylabel('Phase (\pi rad)');
        grid on;

        % Find the max magnitude for the current band
        maxMagnitudes(i) = max(abs(Y));
    end
    
    % Determine the octave for the region based on the band with the highest peak
    [~, maxBandIdx] = max(maxMagnitudes); % Index of the band with the highest peak
    regionOctave(regionIdx) = maxBandIdx + 1; % Octave is band index + 1
    
    % Display the result
    disp(['Region ', num2str(regionIdx), ' belongs to Octave ', num2str(regionOctave(regionIdx))]);
end

validationPassed = true; % Flag to track validation status

for regionIdx = 1:3
    fprintf('Validating Region %d...\n', regionIdx);
    fprintf('Computed Octave: %d | Expected Octave: %d\n', regionOctave(regionIdx), expectedOctaves(regionIdx));
    
    if regionOctave(regionIdx) == expectedOctaves(regionIdx)
        fprintf('Region %d: Validation Passed.\n', regionIdx);
    else
        fprintf('Region %d: Validation Failed!\n', regionIdx);
        validationPassed = false;
    end
end

if validationPassed
    disp('All regions passed validation!');
else
    disp('Some regions failed validation. Check the filtering and analysis.');
end

%% 5.3e)
% Define filter length and sampling frequency
L = length(h); % Filter length
transientDuration = L / fs; % Transient duration in seconds

fprintf('Transient Duration: %.4f seconds\n', transientDuration);

% Analyze transient effects for each filter in each region
for regionIdx = 1:3
    currentRegion = regions{regionIdx}; % Select the current region
    filterOutputs = zeros(numBands, length(currentRegion));
    
    % Apply each filter to the current region
    for i = 1:numBands
        % Calculate center frequency and normalized cutoff
        fc = (bands{i, "StartingFreq_Hz_"} + bands{i, "EndingFreq_Hz_"}) / 2;
        wc = 2 * pi * fc / fs; % Convert to normalized frequency

        % Call HammingNorm to generate frequency response
        H = HammingNorm(wc, L, N); % Obtain frequency response
        h = ifft(ifftshift(H), 'symmetric'); % Filter coefficients
        h = h(1:L); % Trim to filter length

        % Filter the signal
        filterOutputs(i, :) = filter(h, 1, currentRegion);

        % Plot transient effects
        figure(regionIdx + 20); % Separate figure for each region
        subplot(numBands, 1, i);
        plot(((0:length(currentRegion)-1) / fs), filterOutputs(i, :));
        hold on;
        xline(transientDuration, 'r--', 'LineWidth', 1.5); % Mark transient duration
        title(['Filter Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
            '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
        xlabel('Time (s)');
        xlim([0 0.05])
        ylabel('Amplitude');
        sgtitle(['First 5 Secconds of Region ', num2str(regionIdx)]);
        grid on;
        hold off;
    end
end

% The transient duration lasts roughly 12.6 ms which can be seen in the
% plots indicated by the red vertical line. The duration is the same for
% each filter in the bank.
%%  Project Functions

% Simple Band Pass Filter (4.1a)

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


% Width of pass band (4.1b)

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

% Better Bandpass filter

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