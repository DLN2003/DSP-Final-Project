% Authored By: Giovanni Gutierrez
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 4.2 A Better BPF
%
% This code runs a similar analysis to 4.1c), providing
% context for the usefulness of filter length in the design
% process. The passband is wider than in the simple BPF 
% design which is a cost of reducing the noise level with 
% our Hamming window.
%


%% 4.2b) Find the passband width at a threshold level of 50%. Then plot the graphs for two more BPFs at L = 21 and L = 81 and get their passbands.

% 50% of the maximum 10.88 is 5.44, so this defines our threshold. 

% This threshold can be calculated by using the max function
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
