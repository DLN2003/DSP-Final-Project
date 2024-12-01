%% 4.2d) Use frequency response to explain why the filter only passes at cutoff

% Observing the frequency response of the 41-length filter from 4.2a) we
% can see that the response is only significant at the 0.25*pi frequency.
% The response dies down to very shallow side lobes outside of this
% passband. If you observe the output equation from 4.2c) it is easy to see
% that the magnitude of the response for an input centered anywhere in the
% frequency range of the filter will be either amplified of attentuated by
% the filter. If an input component falls within the passband, it will
% be dominant in the output. However, if it falls within the stop-band, it
% will be greatly attenuated and not be greatly represented in the output.
