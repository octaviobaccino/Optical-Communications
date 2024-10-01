function [o_data_s] = optical_RX(i_rx_s, i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    
    config.M = 16;
    config.BR = 64e9;
    config.pulse_shaping_ntaps = 61;
    config.rolloff = 0.1;
    
    config.down_phase = 0;
    
    config.polarization_swap = false;
    config.polarity_swap_hi = false;
    config.polarity_swap_hq = false;
    config.polarity_swap_vi = false;
    config.polarity_swap_vq = false;
    config.rotation_h = 0*pi/2;
    config.rotation_v = 0*pi/2;  
    
    config.enable_plots_rx = 0;
    
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
    
    M = config.M;
    BR = config.BR;
    down_phase = config.down_phase;
    pulse_shaping_ntaps = config.pulse_shaping_ntaps;
    rolloff = config.rolloff;
    
    polarization_swap = config.polarization_swap;
    polarity_swap_hi = config.polarity_swap_hi;
    polarity_swap_hq = config.polarity_swap_hq;
    polarity_swap_vi = config.polarity_swap_vi;
    polarity_swap_vq = config.polarity_swap_vq;
    rotation_h = config.rotation_h;
    rotation_v = config.rotation_v;   
    enable_plots_rx = config.enable_plots_rx;
    
    fs_rx = 2*BR;

    rx_h = i_rx_s.ch_out_h;
    rx_v = i_rx_s.ch_out_v;
    
    h_mf = root_raised_cosine(BR/2, fs_rx, rolloff, pulse_shaping_ntaps, 0);
    
    %--------------------------%
    %         PROCESS
    %--------------------------%
    
    % Filtrado con matched filter
    rx_h = filter(h_mf, 1, rx_h);
    rx_v = filter(h_mf, 1, rx_v);
    
    % Downsampling
    rx_h = rx_h(1+down_phase:2:end);
    rx_v = rx_v(1+down_phase:2:end);
    
    %---------------%
    % 4D CORRELATOR
    %---------------%

    rx_h_ref = rx_h;
    rx_v_ref = rx_v;
    
    if polarization_swap
        rx_h_aux = rx_h;
        rx_v_aux = rx_v;
        rx_h = rx_v_aux;
        rx_v = rx_h_aux;
    end
    
    if polarity_swap_hi
        rx_h = -real(rx_h) + 1j*imag(rx_h);
    end
    
    if polarity_swap_hq
        rx_h = real(rx_h) - 1j*imag(rx_h);
    end
    
    if polarity_swap_vi
        rx_v = -real(rx_v) + 1j*imag(rx_v);
    end
    
    if polarity_swap_vq
        rx_v = real(rx_v) - 1j*imag(rx_v);
    end
    
    if rotation_h ~= 0
        rx_h = rx_h.*exp(1j*rotation_h);
    end
    
    if rotation_v ~= 0
        rx_v = rx_v.*exp(1j*rotation_v);
    end    
    
    %% CORRELADOR 4D
    
    LCORR=1000;
    xh_trim = rx_h_ref(end-LCORR+1:end);
    xv_trim = rx_v_ref(end-LCORR+1:end);
    yh_trim = rx_h(end-LCORR+1:end);
    yv_trim = rx_v(end-LCORR+1:end);

    corr_hi_hi = xcorr(real(xh_trim), real(yh_trim));
    corr_hi_hq = xcorr(real(xh_trim), imag(yh_trim));
    corr_hi_vi = xcorr(real(xh_trim), real(yv_trim));
    corr_hi_vq = xcorr(real(xh_trim), imag(yv_trim));
    
    corr_hq_hi = xcorr(imag(xh_trim), real(yh_trim));
    corr_hq_hq = xcorr(imag(xh_trim), imag(yh_trim));
    corr_hq_vi = xcorr(imag(xh_trim), real(yv_trim));
    corr_hq_vq = xcorr(imag(xh_trim), imag(yv_trim));
    
    corr_vi_hi = xcorr(real(xv_trim), real(yh_trim));
    corr_vi_hq = xcorr(real(xv_trim), imag(yh_trim));
    corr_vi_vi = xcorr(real(xv_trim), real(yv_trim));
    corr_vi_vq = xcorr(real(xv_trim), imag(yv_trim));
    
    corr_vq_hi = xcorr(imag(xv_trim), real(yh_trim));
    corr_vq_hq = xcorr(imag(xv_trim), imag(yh_trim));
    corr_vq_vi = xcorr(imag(xv_trim), real(yv_trim));
    corr_vq_vq = xcorr(imag(xv_trim), imag(yv_trim));
    
    if enable_plots_rx
        figure
        subplot 441
        plot(corr_hi_hi);
        grid on
        subplot 442
        plot(corr_hi_hq);
        grid on
        subplot 443
        plot(corr_hi_vi);
        grid on
        subplot 444
        plot(corr_hi_vq);
        grid on    
        subplot 445
        plot(corr_hq_hi);
        grid on
        subplot 446
        plot(corr_hq_hq);
        grid on
        subplot 447
        plot(corr_hq_vi);
        grid on
        subplot 448
        plot(corr_hq_vq);
        grid on
        subplot 449
        plot(corr_vi_hi);
        grid on
        subplot (4,4,10)
        plot(corr_vi_hq);
        grid on
        subplot (4,4,11)
        plot(corr_vi_vi);
        grid on
        subplot (4,4,12)
        plot(corr_vi_vq);
        grid on
        subplot (4,4,13)
        plot(corr_vq_hi);
        grid on
        subplot (4,4,14)
        plot(corr_vq_hq);
        grid on
        subplot (4,4,15)
        plot(corr_vq_vi);
        grid on    
        subplot (4,4,16)
        plot(corr_vq_vq);
        grid on
    end
    
    % Matriz de correlaciones 4x4
        
    M_corr = [corr_hi_hi corr_hi_hq corr_hi_vi corr_hi_vq...
        corr_hq_hi corr_hq_hq corr_hq_vi corr_hq_vq...
        corr_vi_hi corr_vi_hq corr_vi_vi corr_vi_vq...
        corr_vq_hi corr_vq_hq corr_vq_vi corr_vq_vq];

    M_aux = zeros(1,16);
    
    for idx1 = 1:16
        M_aux(idx1) = max(M_corr(:,idx1));
        if M_aux(idx1) < 2000
            M_aux(idx1) = min(M_corr(:,idx1));
            if M_aux(idx1) > -2000
                M_aux(idx1) = 0;
            end
        end
        if M_aux(idx1) ~= 0
            M_aux(idx1) = M_aux(idx1)/abs(M_aux(idx1));
        end
    end
    
    M_detection = [M_aux(1) M_aux(2) M_aux(3) M_aux(4);...
        M_aux(5) M_aux(6) M_aux(7) M_aux(8);...
        M_aux(9) M_aux(10) M_aux(11) M_aux(12);...
        M_aux(13) M_aux(14) M_aux(15) M_aux(16)];

    if M_detection(1,3) ~= 0 || M_detection(1,4) ~= 0 % detect_polarization_swap
        rx_h_aux = rx_h;
        rx_v_aux = rx_v;
        rx_h = rx_v_aux;
        rx_v = rx_h_aux;
    end
    
    if M_detection(1,2) ~= 0 || M_detection(1,4) ~= 0 % detect_rotation_h
        rx_h = imag(rx_h) + 1j*real(rx_h);
    end
    
    if M_detection(3,2) ~= 0 || M_detection(3,4) ~= 0 % detect_rotation_v
        rx_v = imag(rx_v) + 1j*real(rx_v);
    end
    
    if M_detection(1,1) < 0 || M_detection(1,2) < 0 || M_detection(1,3) < 0 || M_detection(1,4) < 0 % detect_polarity_swap_hi
        rx_h = -real(rx_h) + 1j*imag(rx_h);
    end
    
    if M_detection(2,1) < 0 || M_detection(2,2) < 0 || M_detection(2,3) < 0 || M_detection(2,4) < 0 % detect_polarity_swap_hq
        rx_h = real(rx_h) - 1j*imag(rx_h);
    end    

    if M_detection(3,1) < 0 || M_detection(3,2) < 0 || M_detection(3,3) < 0 || M_detection(3,4) < 0 % detect_polarity_swap_vi
        rx_v = -real(rx_v) + 1j*imag(rx_v);
    end

    if M_detection(4,1) < 0 || M_detection(4,2) < 0 || M_detection(4,3) < 0 || M_detection(4,4) < 0 % detect_polarity_swap_vq
        rx_v = real(rx_v) - 1j*imag(rx_v);
    end    
    
    if enable_plots_rx
        xh_trim = rx_h_ref(end-LCORR+1:end);
        xv_trim = rx_v_ref(end-LCORR+1:end);
        yh_trim = rx_h(end-LCORR+1:end);
        yv_trim = rx_v(end-LCORR+1:end);

        corr_hi_hi = xcorr(real(xh_trim), real(yh_trim));
        corr_hi_hq = xcorr(real(xh_trim), imag(yh_trim));
        corr_hi_vi = xcorr(real(xh_trim), real(yv_trim));
        corr_hi_vq = xcorr(real(xh_trim), imag(yv_trim));

        corr_hq_hi = xcorr(imag(xh_trim), real(yh_trim));
        corr_hq_hq = xcorr(imag(xh_trim), imag(yh_trim));
        corr_hq_vi = xcorr(imag(xh_trim), real(yv_trim));
        corr_hq_vq = xcorr(imag(xh_trim), imag(yv_trim));

        corr_vi_hi = xcorr(real(xv_trim), real(yh_trim));
        corr_vi_hq = xcorr(real(xv_trim), imag(yh_trim));
        corr_vi_vi = xcorr(real(xv_trim), real(yv_trim));
        corr_vi_vq = xcorr(real(xv_trim), imag(yv_trim));

        corr_vq_hi = xcorr(imag(xv_trim), real(yh_trim));
        corr_vq_hq = xcorr(imag(xv_trim), imag(yh_trim));
        corr_vq_vi = xcorr(imag(xv_trim), real(yv_trim));
        corr_vq_vq = xcorr(imag(xv_trim), imag(yv_trim));

        figure
        subplot 441
        plot(corr_hi_hi);
        grid on
        subplot 442
        plot(corr_hi_hq);
        grid on
        subplot 443
        plot(corr_hi_vi);
        grid on
        subplot 444
        plot(corr_hi_vq);
        grid on    
        subplot 445
        plot(corr_hq_hi);
        grid on
        subplot 446
        plot(corr_hq_hq);
        grid on
        subplot 447
        plot(corr_hq_vi);
        grid on
        subplot 448
        plot(corr_hq_vq);
        grid on
        subplot 449
        plot(corr_vi_hi);
        grid on
        subplot (4,4,10)
        plot(corr_vi_hq);
        grid on
        subplot (4,4,11)
        plot(corr_vi_vi);
        grid on
        subplot (4,4,12)
        plot(corr_vi_vq);
        grid on
        subplot (4,4,13)
        plot(corr_vq_hi);
        grid on
        subplot (4,4,14)
        plot(corr_vq_hq);
        grid on
        subplot (4,4,15)
        plot(corr_vq_vi);
        grid on    
        subplot (4,4,16)
        plot(corr_vq_vq);
        grid on
    end
    
    % Slicer
    ak_hat_h = my_slicer(rx_h, M);
    ak_hat_v = my_slicer(rx_v, M);
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%
    
    o_data_s.rx_h = rx_h;
    o_data_s.rx_v = rx_v;
    o_data_s.ak_hat_h = ak_hat_h;
    o_data_s.ak_hat_v = ak_hat_v;
        
end

