%% Gain Imbalance sweep

clear; close all

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% General configuration

% -- General --

config.BR = 64e9;                               % Baudrate
config.M = 16;                                  % Modulacion (QPSK/QAM-16/QAM-64...)

% -- TX --

config.config_tx.Lsymbs = 800e3;                % Cantidad de simbolos transmitidos
config.config_tx.rolloff = 0.1;                 % Exceso de ancho de banda del pulse shaping
config.config_tx.pulse_shaping_ntaps = 51;      % Taps del pulse shaping
config.config_tx.pulse_shaping_type = 0;

config.config_tx.skew = 0;
config.config_tx.p0_dbm = 10;                        
config.config_tx.lw = 0; 
config.config_tx.delta_f = 0;    
config.config_tx.Vpi = 6;
config.config_tx.swing = 0.5;               % max = 2
config.config_tx.ER_inner_db = inf;
config.config_tx.ER_outer_db = inf;         % min = 20
% config.config_tx.phase_error = 0/180*pi;

% -- CH --

config.config_ch.N = 2;                         % Tasa de sobremuestreo del canal
% config.config_ch.osnr_db = 20;                % optical SNR
config.config_ch.pulse_shaping_ntaps = 81;      % Taps del filtro del TX
config.config_ch.rolloff = 0.1;                 % Rolloff del TX
config.config_ch.f_sop = 0*50e3;                  % Frecuencia de rotacion del SOP

% -- RX -- 

config.config_rx.down_phase = 0;                % Fase de downsampling
config.config_rx.tap_leak = 1e-5;               % tap leakeage
config.config_rx.target_agc = 0.3;              % ganancia objetivo del AGC
config.config_rx.eq_taps = 31;                  % numero de taps del ecualizador MIMO
config.config_rx.step_cma = 2^-9;               % step del CMA
config.config_rx.step_dd = 2^-9;                % step del DD


config.config_rx.polarization_swap = false;
config.config_rx.polarity_swap_hi = false;
config.config_rx.polarity_swap_hq = false;
config.config_rx.polarity_swap_vi = false;
config.config_rx.polarity_swap_vq = false;
config.config_rx.rotation_h = 0*pi/2;
config.config_rx.rotation_v = 0*pi/2;  
config.config_rx.enable_plots_rx = 0;

%% Sweep

phase_v = [0 2 4 6 8 10 12 14];
n_phase = length(phase_v);

ber_max = 2e-2;
ber_min = 1e-3;
n_ber = 10;
linespace_v = linspace(log10(ber_min), log10(ber_max), n_ber);
theo_ber_v = 10.^(linespace_v);

% Save configuration
file = [out_dir, 'cfg.mat'];
save(file,'phase_v','theo_ber_v','config');

out_c = cell(n_ber, 1); 

% Loop

M = config.M;
BR = config.BR;

for idx_phase = 1:n_phase
    
    phase_error = phase_v(idx_phase);

    osnr_db_v = get_osnr_from_theo_ber(theo_ber_v,M,BR);
    
    config.config_tx.phase_error = phase_error*pi/180;

    for idx_ber = 1:n_ber

        osnr_db = osnr_db_v(idx_ber);
        config_aux = config;
        config_aux.config_ch.osnr_db = osnr_db;

        fprintf('- Running phase-%d/%d-ber-%d/%d ...\n', ...
            idx_phase, n_phase, idx_ber, n_ber)

        out_c{idx_ber} = optical_simulator_v3(config_aux);

    end

    name = sprintf('phase_%d',phase_error);
    file = [out_dir, 'out_',name,'.mat'];
    save(file, 'out_c');
    
end

