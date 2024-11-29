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
