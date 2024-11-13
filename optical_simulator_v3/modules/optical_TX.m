function [o_data_s] = optical_TX(i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    config.BR = 64e9;                   % Baudrate
    config.M = 16;                      % Cantidad de niveles de la modulacion
    config.Lsymbs = 100e3;              % Cantidad de simbolos transmitidos
    config.rolloff = 0.1;               % Rolloff del filtro conformador
    config.pulse_shaping_ntaps = 61;    % Taps del filtro conformador
    config.pulse_shaping_type = 0;      % 0: RRC, 1: RC
    config.skew = 0;                    % skew como fraccion del periodo de simbolo
    config.p0_dbm = 20;                 % potencia en dBm
    config.lw = 0;                      % linewidth
    config.delta_f = 0;                 % frequency offset
    config.Vpi = 6;
    config.swing = 0.8; % Se define como el voltaje pico a pico a la entrada de los electrodos del MZI como porcentaje del Vpi. Max swing=2
    config.ER_inner_db = inf;
    config.ER_outer_db = inf;
    config.phase_error = 0/180*pi;
    
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

    %%
    %--------------------------%
    %         VARIABLES
    %--------------------------%
    
    BR = config.BR;
    M = config.M;
    Lsymbs = config.Lsymbs;
    rolloff = config.rolloff;
    pulse_shaping_ntaps = config.pulse_shaping_ntaps;
    pulse_shaping_type = config.pulse_shaping_type;
    skew = config.skew;
    p0_dbm = config.p0_dbm;
    lw = config.lw;
    delta_f = config.delta_f;
    
    % Transmisor siempre a tasa 2
    fs_tx = BR*2;

    %%
    %--------------------------%
    %         PROCESS
    %--------------------------%

    % QAM Symbols generation 
    dec_labels_v = randi([0 M-1], Lsymbs, 1);
    dec_labels_h = randi([0 M-1], Lsymbs, 1);

    tx_symbs_v = qammod(dec_labels_v,M);
    tx_symbs_h = qammod(dec_labels_h,M);

    % Upsampling to change sampling rate
    xup_v = 2*upsample(tx_symbs_v, 2);
    xup_h = 2*upsample(tx_symbs_h, 2);

    % Pulse shaping con RRC-RC
    if pulse_shaping_type==0
        h_ps = root_raised_cosine(BR/2, fs_tx, rolloff, pulse_shaping_ntaps, 0);
    else
        h_ps = raised_cosine(BR/2, fs_tx, rolloff, pulse_shaping_ntaps, 0);
    end
 
    % Signal shaping
    yup_v = filter(h_ps,1,xup_v);
    yup_h = filter(h_ps,1,xup_h);
    
    clear xup_v
    clear xup_h
    
    %%
    %--------------------------%
    %           TOSA
    %--------------------------%
    
    % DAC (amp. normalizer)
    
    Vp = 0.5; % Maxima excursion del DAC
    
    % Normalizacion (mal hecha, corregir en algun momento)
    yup_h = Vp*real(yup_h)/var(real(yup_h)) + Vp*1j*imag(yup_h)/var(real(yup_h));
    yup_v = Vp*real(yup_v)/var(real(yup_v)) + Vp*1j*imag(yup_v)/var(real(yup_v));
    
    % Delay
    
    skew_s = skew/BR; % paso el skew a segundos
    
    yup_h = skew_inserter(yup_h, skew_s, fs_tx); % skew en segundos
    yup_v = skew_inserter(yup_v, skew_s, fs_tx);
    
    % MZM - Driver Response (para este simulador solo es una ganancia)
    
    MZM_gain = 5;
    
    h_mzm = 1; % filtro para representar respuesta en frecuencia del MZM
    
    yup_h = MZM_gain*yup_h;
    yup_h = filter(h_mzm,1,yup_h);

    yup_v = MZM_gain*yup_v;
    yup_v = filter(h_mzm,1,yup_v);
    
    % Dual Polarization MZM (con su laser)
    
    y_mzi_h = MZM_function(yup_h,config);
    y_mzi_v = MZM_function(yup_v,config);
    
%     laser_generator(potencia,LW,frec_offset,fs,N,Lsymbs)

%     p_laser = laser_generator(p0,lw,delta_f,fs_tx,2,Lsymbs); % p_laser =  
    
    %%
    %--------------------------%
    %         OUTPUT
    %--------------------------%

    o_data_s.tx_out_v = y_mzi_v;
    o_data_s.tx_out_h = y_mzi_h;
    
    o_data_s.ak_h = tx_symbs_h;
    o_data_s.ak_v = tx_symbs_v;
    
    o_data_s.h_ps = h_ps;
    
end