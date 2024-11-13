
function [ipr] = my_rcosine(BR, fs, rolloff, NTAPS, t0)
%BR = 33e9; fs=33e9; rolloff=0.20001; NTAPS=256; phase=1/33e9/2;
rolloff = rolloff + 0.0001;
Ts = 1/fs;
T = 1/BR;

t = [-floor(NTAPS/2):1:floor(NTAPS/2)].*Ts + t0;
t_norm = t./T;
ipr = sinc(t_norm).*(cos(pi*rolloff.*t_norm))./(1-(2*rolloff.*t_norm).^2);

end