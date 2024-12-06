%% 5.3a) (by David Needens Nov. 2024)
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
