% Authored By: David Needens
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 5.3 Piano Octave Decoding
% 
% This MATLAB code plots the time-domain outputs of multiple bandpass 
% filters, each corresponding to a specific frequency range, in a subplot 
% format, with one plot per filter displaying its respective filtered 
% signal.
%

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
