function [o_data_s] = optical_CH(i_ch_s, i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    
    config.BR = 64e9;         % Cantidad de niveles de la modulacion
    config.M = 16;            % Orden de modulacion
    config.N = 2;             % Factor de sobremuesteo
    
    config.osnr_db = 20;      % OSNR [dB]
    
    config.pulse_shaping_ntaps = 81; % Parametros del filtro TX/RX
    config.rolloff = 0.1;
    
    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%
    
    fn = fieldnames(i_config_s);
    for k = 1:numel(fn)
        if isfield(config,(fn{k}))==1
            config.(fn{k})= i_config_s.(fn{k});
        else
            error("%s: Parametro del simulador no valido", fn{k})
        end
    end
    
    %--------------------------%
    %         VARIABLES
    %--------------------------%
    
    osnr_db = config.osnr_db;           % Eb/No en dB
    M = config.M;                       % Orden de modulacion
    N = config.N;                       % Factor de sobremuesteo
    pulse_shaping_ntaps = config.pulse_shaping_ntaps;
    rolloff = config.rolloff;
    BR = config.BR;
    tx_out_v = i_ch_s.tx_out_v;         % Senial de entrada al canal V
    tx_out_h = i_ch_s.tx_out_h;         % Senial de entrada al canal H
    
    % Frecuencia de muestreo del canal (N veces la tasa 2 del transmisor)
    fs_ch = 2*BR*N;
    
    %--------------------------%
    %         PROCESS
    %--------------------------%
    
    % Upsampling del canal
    xup_h = N*upsample(tx_out_h,N);
    xup_v = N*upsample(tx_out_v,N);
    
    % Filtro para simular ancho de banda del transmisor (frecuencia de corte 1.5 veces la frecuencia de nyquist de la senial transmitida)
    h_tx = root_raised_cosine(BR/2*1.5, fs_ch, rolloff, pulse_shaping_ntaps, 0);
    
    yup_h = filter(h_tx, 1, xup_h); 
    yup_v = filter(h_tx, 1, xup_v);
   
    % Generacion de ruido
    snr_db = osnr_db - 10*log10(BR/12.5e9);
    snr = 10^(snr_db/10);
    Ps = var(yup_h(10e3:50e3))*2*N; % Por que se multiplica por 2N?
    Pn = Ps/snr;
    
    n_h = sqrt(Pn/2).*(randn(size(yup_h)) + 1j.*randn(size(yup_h)));
    n_v = sqrt(Pn/2).*(randn(size(yup_v)) + 1j.*randn(size(yup_v)));
    
    % Agrego ruido a la senial
    ych_h = yup_h + n_h;
    ych_v = yup_v + n_v;

    % Filtro para simular ancho de banda del receptor
    h_rx = root_raised_cosine(BR/2*1.5, fs_ch, rolloff, pulse_shaping_ntaps, 0);
    
    yrx_h = filter(h_rx,1,ych_h);
    yrx_v = filter(h_rx,1,ych_v);
    
    % Downsampling del canal para entrar al receptor
    yrx_h = yrx_h(1:N:end);
    yrx_v = yrx_v(1:N:end);
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%
    
    % Salida pre filtro TX
    o_data_s.xup_ch_h = xup_h;
    o_data_s.xup_ch_v = xup_v;
    
    % Salida post filtro TX
    o_data_s.yup_ch_h = yup_h;
    o_data_s.yup_ch_v = yup_v;
    
    % Salida con ruido
    o_data_s.ych_h = ych_h;
    o_data_s.ych_v = ych_v;
    
    % Salida al receptor
    
    o_data_s.ch_out_h = yrx_h;
    o_data_s.ch_out_v = yrx_v;
    
end