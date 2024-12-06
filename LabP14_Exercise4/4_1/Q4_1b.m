
% Authored By: Giovanni Gutierrez
%
% Last Modified: Nov. 2024
%
% Reference 4 Lab Exercises: 4.1 Bandpass Filter Design
%
% Knowing the range of frequencies that the BPF will 
% pass represent the accepted frequencies of the filter.
% With a chosen threshold this code utilizes a project fxn
% designed to measure the width of the passband.
%


%% 4.1b) Using the 0.5 level passband, find width of passband 

% Use Project function to define the pass band width
Passband = PBWidth(H(idx:end),w(idx:end),0.5);

% Display passband width
fprintf('Passband width @ the 0.5 level is approximately %.4f radians for the L = 40 filter\n', Passband);

% We can imagine this width drawn on our first plot, in 4.1a), between the sides of
% the main lobe that is centered at wc. It can give us a range of
% frequencies accepted by our passband using the cutoff frequency where
% Range of accepted w = wc +- PassBand_Width/2 
