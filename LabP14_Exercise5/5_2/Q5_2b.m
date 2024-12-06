%% 5.2b)
wc = BP_Filters.CenterFreq_omegaC_;
N = 2048; % # of Points for the DFT
w = -pi: 2*pi/N: pi - 2*pi/N;%Frequency range
% Filters for octaves 1-7. L was determined through trial and error and was
% chosen based on which L generated a bandwidth closest to the design
% bandwidth for each octave.
% L = [519, 251, 128,63,32,17,9];% Values originally found for 0.5 stopband
L = [943, 546, 246,133, 67, 34, 17]; % Values with a 0.01 stopband
