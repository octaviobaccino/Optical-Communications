%-------------------------------------------------------------------------%
%   FUNDACION FULGOR
%-------------------------------------------------------------------------%
%   OPTICAL SIMULATOR (version 1)
%   Autor: Octavio Baccino
%-------------------------------------------------------------------------%

function odata = optical_simulator_v1(i_config_s)

close all;

%------------------%
% DEFAULT SETTINGS
%------------------%

% General

config.BR = 64e9;
config.M = 16;

% TX

config.config_tx.BR = config.BR;
config.config_tx.M = config.M;
config.config_tx.Lsymbs = 500e3;
config.config_tx.rolloff = 0.1;
config.config_tx.pulse_shaping_ntaps = 61;
config.config_tx.pulse_shaping_type = 0;

% CH

config.config_ch.BR = config.BR;
config.config_ch.M = config.M;
config.config_ch.N = 2;
config.config_ch.osnr_db = 22;
config.config_ch.pulse_shaping_ntaps = 81;
config.config_ch.rolloff = 0.1;

% RX 

config.config_rx.BR = config.BR;
config.config_rx.M = config.M;
config.config_rx.pulse_shaping_ntaps = config.config_tx.pulse_shaping_ntaps;
config.config_rx.rolloff = config.config_tx.rolloff;
config.config_rx.down_phase = 0;

config.config_rx.polarization_swap = false;
config.config_rx.polarity_swap_hi = false;
config.config_rx.polarity_swap_hq = false;
config.config_rx.polarity_swap_vi = false;
config.config_rx.polarity_swap_vq = false;
config.config_rx.rotation_h = 0*pi/2;
config.config_rx.rotation_v = 0*pi/2;  
config.config_rx.enable_plots_rx = 0;

%---------------------%
% OVERWRITE PARAMETERS
%---------------------%

if nargin > 0
    config = overwrite_parameters(i_config_s, config);
end

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

odata = rx_odata;

%-----------%
% BER CHECK
%-----------%

snr_db = config.config_ch.osnr_db - 10*log10(config.BR/12.5e9);
ebno = snr_db-10*log10(log2(config.M));

ber_theo = berawgn(ebno, 'QAM', config.M);

ak_h = tx_odata.ak_h;
ak_v = tx_odata.ak_v;

ak_hat_h = rx_odata.ak_hat_h;
ak_hat_v = rx_odata.ak_hat_v;

[ber_est_h, n_errors_h] = my_ber_checker(ak_hat_h,ak_h,config.M,'auto');
[ber_est_v, n_errors_v] = my_ber_checker(ak_hat_v,ak_v,config.M,'auto');
ber_est = (ber_est_h + ber_est_v)/2;

odata.ber_theo = ber_theo;
odata.ber_est = ber_est;

end



