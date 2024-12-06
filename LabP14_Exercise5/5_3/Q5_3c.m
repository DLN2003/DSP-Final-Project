%% 5.3c) (by David Needens Nov. 2024)
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
