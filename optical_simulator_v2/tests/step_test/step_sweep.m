%% Sweep step

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

% config.BR = 64e9;                               % Baudrate
config.M = 16;                                  % Modulacion (QPSK/QAM-16/QAM-64...)

% -- TX --

config.config_tx.Lsymbs = 300e3;                % Cantidad de simbolos transmitidos
config.config_tx.rolloff = 0.1;                 % Exceso de ancho de banda del pulse shaping
config.config_tx.pulse_shaping_ntaps = 81;      % Taps del pulse shaping
config.config_tx.pulse_shaping_type = 0;

% -- CH --

config.config_ch.N = 2;                         % Tasa de sobremuestreo del canal
% config.config_ch.osnr_db = 20;                % optical SNR
config.config_ch.pulse_shaping_ntaps = 81;      % Taps del filtro del TX
config.config_ch.rolloff = 0.1;                 % Rolloff del TX
config.config_ch.f_sop = 0*50e3;                  % Frecuencia de rotacion del SOP

% -- RX -- 

config.config_rx.down_phase = 0;                % Fase de downsampling
config.config_rx.tap_leak = 1e-2;               % tap leakeage
config.config_rx.target_agc = 0.3;              % ganancia objetivo del AGC
config.config_rx.eq_taps = 31;                  % numero de taps del ecualizador MIMO
config.config_rx.step_cma = 2e-4;               % step del CMA
% config.config_rx.step_dd = 2e-4;                % step del DD


config.config_rx.polarization_swap = false;
config.config_rx.polarity_swap_hi = false;
config.config_rx.polarity_swap_hq = false;
config.config_rx.polarity_swap_vi = false;
config.config_rx.polarity_swap_vq = false;
config.config_rx.rotation_h = 0*pi/2;
config.config_rx.rotation_v = 0*pi/2;  
config.config_rx.enable_plots_rx = 0;

%% Sweep

BR_v = [32e9 64e9];
n_BR = length(BR_v);

step_dd_v = [2e-5 2e-6 2e-7 2e-8];
n_steps = length(step_dd_v);

ber_max = 1e-1;
ber_min = 5e-4;
n_ber = 10;
linespace_v = linspace(log10(ber_min), log10(ber_max), n_ber);
theo_ber_v = 10.^(linespace_v);

% Save configuration
file = [out_dir, 'cfg.mat'];
save(file,'BR_v','step_dd_v','theo_ber_v','config');

out_c = cell(n_ber, 1); 

% Loop

M = config.M;

for idx_BR = 1:n_BR
    
    BR = BR_v(idx_BR);
    
    out_dir_BR = [out_dir(1:end-1), '\'];
    name = sprintf('out_BR_%2d/', BR/1e9);
    out_dir_BR = [out_dir_BR, name];

    if ~exist(out_dir_BR,'dir')
        mkdir(out_dir_BR);
    end
    
    osnr_db_v = get_osnr_from_theo_ber(theo_ber_v,M,BR);
    
    config.BR = BR;
    
    for idx_step = 1:n_steps
        
        step_dd = step_dd_v(idx_step);
        
        config.config_rx.step_dd = step_dd;
        
        for idx_ber = 1:n_ber

            config.config_ch.osnr_db = osnr_db_v(idx_ber);

            fprintf('- Running BR-%d/%d-step-%d/%d-ber-%d/%d ...\n', ...
                idx_BR,n_BR, idx_step, n_steps, idx_ber, n_ber)

            out_c{idx_ber} = optical_simulator_v2(config);

        end
        
        name = sprintf('step_dd_%e',step_dd);
        file = [out_dir_BR, 'out_',name,'.mat'];
        save(file, 'out_c');
        
    end
    
end
