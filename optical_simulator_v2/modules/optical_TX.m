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
    
    BR = config.BR;
    M = config.M;
    Lsymbs = config.Lsymbs;
    rolloff = config.rolloff;
    pulse_shaping_ntaps = config.pulse_shaping_ntaps;
    pulse_shaping_type = config.pulse_shaping_type;
    
    % Transmisor siempre a tasa 2
    fs_tx = BR*2;

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
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%

    o_data_s.tx_out_v = yup_v;
    o_data_s.tx_out_h = yup_h;
    
    o_data_s.ak_h = tx_symbs_h;
    o_data_s.ak_v = tx_symbs_v;
    
    o_data_s.h_ps = h_ps;
    
end