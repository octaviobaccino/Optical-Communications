%-------------------------------------------------------------------------%
% Filename		: root_raised_cosine.m
% Programmer    : Patricio Reus Merlo
% Created on	: 08/01/2023
% Description 	: Raised cosine filter
%-------------------------------------------------------------------------%

function [h_v, n_taps] = root_raised_cosine(fc, fs, rolloff, n_taps, t0)
    
    if nargin < 5
        t0 = 0;
    end
    
    Ts = 1/fs;
    rolloff = rolloff + 0.0001;
    T  = 1 / fc ;
    
    % Force to odd
    n_taps = round_odd(n_taps);
    
    % Time vector
    t_v = ((- (n_taps - 1) / 2 : 1 : (n_taps - 1) / 2) .* Ts + t0);
    tn_v = t_v * 2/T;
    
    % Filter taps
    h_v = (sin(pi .* (1 - rolloff) .* tn_v) ...
                + 4 * rolloff .* tn_v .* cos(pi*(1 + rolloff) .* tn_v)) ...
                            ./ (pi .* tn_v .* (1 - (4 * rolloff .* tn_v).^2));

    h_v(1 + (n_taps - 1) / 2) = (1 + rolloff .* (4 ./ pi - 1));
    h_v = h_v ./ sum(h_v);
    
end