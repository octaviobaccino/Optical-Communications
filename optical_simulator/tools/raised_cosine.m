%-------------------------------------------------------------------------%
% Filename		: raised_cosine.m
% Programmer    : Patricio Reus Merlo
% Created on	: 08/01/2023
% Description 	: Raised cosine filter
%-------------------------------------------------------------------------%

function [h_v, n_taps] = raised_cosine(fc, fs, rolloff, n_taps, t0)
    
    if nargin < 5
        t0 = 0;
    end

    rolloff = rolloff + 0.0001; 
    Ts = 1 / fs ;
    T  = 1 / fc ;
    
    % Force to odd
    n_taps = round_odd(n_taps);
    
    % Time vector
    t_v = ((- (n_taps - 1) / 2 : 1 : (n_taps - 1) / 2) .* Ts + t0);
    tn_v = t_v * 2/T;
    
    % Filter taps
    h_v = sinc(tn_v).*(cos(pi .* rolloff .* tn_v )) ...
                                        ./ (1 - (2 * rolloff .* tn_v) .^ 2); 
    h_v = h_v ./ sum(h_v);
end