% Authored By: Cade Boynton
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 5.2 Piano Octave Decoding
%
% 
% This MATLAB code generates the frequency response for each octaves
% Normalized Hanning filter. Here we also display the frequency response
% and center frequency for each filter.
%
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
