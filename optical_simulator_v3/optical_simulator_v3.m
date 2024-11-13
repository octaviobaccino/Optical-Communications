%-------------------------------------------------------------------------%
%   FUNDACION FULGOR
%-------------------------------------------------------------------------%
%   OPTICAL SIMULATOR (version 1)
%   Autor: Octavio Baccino
%-------------------------------------------------------------------------%

function odata = optical_simulator_v3(i_config_s)

close all;

%------------------%
% DEFAULT SETTINGS
%------------------%

% General

config.BR = 64e9;
config.M = 16;

% TX

config.config_tx.Lsymbs = 1000e3;
config.config_tx.rolloff = 0.1;
config.config_tx.pulse_shaping_ntaps = 51;
config.config_tx.pulse_shaping_type = 0;

config.config_tx.skew = 0.15;
config.config_tx.p0_dbm = 10;                        
config.config_tx.lw = 0; 
config.config_tx.delta_f = 0;    
config.config_tx.Vpi = 6;
config.config_tx.swing = 0.5;               % max = 2
config.config_tx.ER_inner_db = inf;
config.config_tx.ER_outer_db = inf; % min = 18
config.config_tx.phase_error = 0/180*pi;

% CH

config.config_ch.N = 2;
config.config_ch.osnr_db = 25;
config.config_ch.pulse_shaping_ntaps = 81;
config.config_ch.rolloff = 0.1;
config.config_ch.f_sop = 0*100e3;

% RX 

config.config_rx.down_phase = 0;
config.config_rx.tap_leak = 1e-5;
config.config_rx.target_agc = 0.3;
config.config_rx.eq_taps = 31;
config.config_rx.step_cma = 2^-9; % Con este step funciona el FSE con CMA (2^-9)
config.config_rx.step_dd = 2^-9;

config.config_rx.polarization_swap = false;
config.config_rx.polarity_swap_hi = false;
config.config_rx.polarity_swap_hq = false;
config.config_rx.polarity_swap_vi = false;
config.config_rx.polarity_swap_vq = false;
config.config_rx.rotation_h = 0*pi/2;
config.config_rx.rotation_v = 0*pi/2;  
config.config_rx.enable_plots_rx = 1;

%---------------------%
% OVERWRITE PARAMETERS
%---------------------%

if nargin > 0
    config = overwrite_parameters(i_config_s, config);
end

% -- Shared Variables -- %

config.config_tx.BR = config.BR;
config.config_tx.M = config.M;
config.config_ch.BR = config.BR;
config.config_ch.M = config.M;
config.config_rx.BR = config.BR;
config.config_rx.M = config.M;
config.config_rx.pulse_shaping_ntaps = config.config_ch.pulse_shaping_ntaps;
config.config_rx.rolloff = config.config_ch.rolloff;

%------------%
% SIMULATION
%------------%

% Transmision de datos

tx_odata = optical_TX(config.config_tx);

ch_odata = optical_CH(tx_odata, config.config_ch);

rx_odata = optical_RX(ch_odata, config.config_rx);

%--------%
% OUTPUT
%--------%

% odata = rx_odata;

%-----------%
% BER CHECK
%-----------%

snr_db = config.config_ch.osnr_db - 10*log10(config.BR/12.5e9);
ebno = snr_db - 10*log10(log2(config.M));

ber_theo = berawgn(ebno, 'QAM', config.M);

ak_h = tx_odata.ak_h;
ak_v = tx_odata.ak_v;

ak_hat_h = rx_odata.ak_hat_h;
ak_hat_v = rx_odata.ak_hat_v;

valid = fix(0.4*length(ak_h)):length(ak_h)-100;

[ber_est_h, n_errors_h] = my_ber_checker(ak_hat_h(valid),ak_h(valid),config.M,'auto');
[ber_est_v, n_errors_v] = my_ber_checker(ak_hat_v(valid),ak_v(valid),config.M,'auto');
ber_est = (ber_est_h + ber_est_v)/2;

odata.ber_theo = ber_theo;
odata.ber_est = ber_est;

end

