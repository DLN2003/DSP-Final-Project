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

maxMagnitude = [218.75, 34.7275, 259.087];

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
        plot(f_band, abs(Y) / maxMagnitude(regionIdx));
        title(['Normalized Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
            '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
        xlabel('Frequency (Hz)');
        xlim([f_band(1) f_band(end)]);
        ylim([0 1]);
        ylabel('Magnitude');
        grid on;

        subplot(numBands, 2, 2*i);
        plot(f_band, angle(Y) / pi);
        title(['Phase of Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
            '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
        xlabel('Frequency (Hz)');
        xlim([f_band(1) f_band(end)]);
        ylim([-1 1]);
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
transientDuration = L / fs * 1000; % Transient duration in milliseconds

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
        figure(regionIdx + 14); % Separate figure for each region
        subplot(numBands, 1, i);
        plot(((0:length(currentRegion)-1) / fs * 1000), filterOutputs(i, :));
        hold on;
        xline(transientDuration, 'r--', 'LineWidth', 1.5); % Mark transient duration
        title(['Filter Output for Band ', num2str(round(bands{i, "StartingFreq_Hz_"})), ...
            '-', num2str(round(bands{i, "EndingFreq_Hz_"})), ' Hz']);
        xlabel('Time (ms)');
        xlim([0 25])
        ylabel('Amplitude');
        sgtitle(['First 25ms of Region ', num2str(regionIdx)]);
        grid on;
        hold off;
    end
end

% The transient duration lasts roughly 12.6 ms which can be seen in the
% plots indicated by the red vertical line. The duration is the same for
% each filter in the bank.
%%  Project Functions


%% Normalized Hanning Bandpass Filter (5.2)

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
