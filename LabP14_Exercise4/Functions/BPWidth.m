% Authored By: Giovanni Gutierrez
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises
%
% This code defines a function that measures the
% width of the passband for a given BPF based on
% the desired threshold.
% 


%% Width of pass band (4.1b)

% Function that finds passband width using the magnitude response, omega,
% and the threshold level for the passband.

function Passband = PBWidth(H,w,th)

% Find frequencies where |H(w)| is above threshold
Hmag = abs(H); % Defintion of magnitude response
H1 = find(Hmag >= th, 1, 'first'); % Find first index where |H(w)| is close to threshold
H2 = find(Hmag >= th, 1, 'last'); % Find last index where |H(w)| is close to threshold

% Convert indices to frequency values
w1 = w(H1); % Frequency at start of passband
w2 = w(H2); % Frequency at end of passband
Passband = w2 - w1; % Width of the passband
end
