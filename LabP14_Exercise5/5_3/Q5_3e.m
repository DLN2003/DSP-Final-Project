% Authored By: David Needens
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 5.3 Piano Octave Decoding
% 
% This code calculates the transient duration for each bandpass filter 
% applied to different signal regions and visualizes the effects by 
% plotting the filter outputs for the first 5 seconds, marking the 
% transient duration with a red dashed line on each plot.
%

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
