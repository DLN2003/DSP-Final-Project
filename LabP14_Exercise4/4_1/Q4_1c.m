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
