function [o_data_s] = optical_RX(i_rx_s, i_config_s)
    
    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%
    
    config.M = 16;
    config.BR = 64e9;
    config.pulse_shaping_ntaps = 61;
    config.rolloff = 0.1;
    config.NRX = 2; % simpre va a ser 2
    config.tap_leak = 1e-2;
    config.target_agc = 0.3;
    config.eq_taps = 31;
    config.step_cma = 2e-3;
    config.step_dd = 2e-3;
    
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
    ak_h = i_rx_s.ak_h; % simbolos transmitidos para calcular RCMA
    ak_v = i_rx_s.ak_v; % simbolos para correaldor
    
    h_mf = root_raised_cosine(BR/2, fs_rx, rolloff, pulse_shaping_ntaps, 0);
    
    %--------------------------%
    %         PROCESS
    %--------------------------%
    
    %%
    %-----%
    % AGC
    %-----%
    
    target_agc = 0.3;
    metric_h = std(rx_h);
    metric_v = std(rx_v);
    rx_agc_h = rx_h/metric_h*target_agc;
    rx_agc_v = rx_v/metric_v*target_agc;
    
    %%
    %---------------%
    % FSE (CMA + DD)
    %---------------%
    
    % Parametros
    eq_taps = config.eq_taps;
    step_cma = config.step_cma;
    step_dd = config.step_dd;
    tap_leak = config.tap_leak;
    Lrx = length(rx_agc_h);
    NRX = 2; % simpre trabaja a tasa 2
    RCMA = sqrt(mean(abs(ak_h(1:10e3)).^4)/mean(abs(ak_h(1:10e3)).^2));
    subsf = 500;
    
    % Ecaulizador MIMO
    
    eq_buffer_h = zeros(eq_taps,1);
    eq_buffer_v = zeros(eq_taps,1);
    heq_00 = zeros(eq_taps,1);
    heq_01 = zeros(eq_taps,1);
    heq_10 = zeros(eq_taps,1);
    heq_11 = zeros(eq_taps,1);
    heq_00((eq_taps+1)/2) = 1;
    heq_11((eq_taps+1)/2) = 1;
    
    % Logeos para plottear
    ak_hat_h = zeros(Lrx/2,1);
    ak_hat_v = zeros(Lrx/2,1);
    yrx_h_log = zeros(Lrx/2,1);
    yrx_v_log = zeros(Lrx/2,1);
    error_h_log = zeros(Lrx/2,1);
    error_v_log = zeros(Lrx/2,1);
    heq_00_log = zeros(Lrx/2/subsf,eq_taps);
    heq_01_log = zeros(Lrx/2/subsf,eq_taps);
    heq_10_log = zeros(Lrx/2/subsf,eq_taps);
    heq_11_log = zeros(Lrx/2/subsf,eq_taps);
    
    for n = 1:Lrx
        
        if n < fix(Lrx/3)
            cma_enable = 1;
        else
            cma_enable = 0;
        end
        
        % Meto nueva muestra en el buffer
        eq_buffer_h(2:end) = eq_buffer_h(1:end-1);
        eq_buffer_h(1) = rx_agc_h(n);
        eq_buffer_v(2:end) = eq_buffer_v(1:end-1);
        eq_buffer_v(1) = rx_agc_v(n);
        
        % Producto MIMO
        yrx_h = sum(eq_buffer_h.*heq_00) + sum(eq_buffer_v.*heq_01);
        yrx_v = sum(eq_buffer_h.*heq_10) + sum(eq_buffer_v.*heq_11);

%         yrx_h_log(n) = yrx_h;
%         yrx_v_log(n) = yrx_v;
        
        if mod(n,NRX)==0
            
            n2 = ceil(n/NRX);
            y_hat_h = my_slicer(yrx_h,M);
            y_hat_v = my_slicer(yrx_v,M);

            ak_hat_h(n2) = y_hat_h;
            ak_hat_v(n2) = y_hat_v;
            
            yrx_h_log(n2) = yrx_h;
            yrx_v_log(n2) = yrx_v;
            
            if cma_enable
                error_h = yrx_h*(abs(yrx_h)-RCMA);
                error_v = yrx_v*(abs(yrx_v)-RCMA);
                step=step_cma;
            else
                error_h = yrx_h - y_hat_h;
                error_v = yrx_v - y_hat_v;
                step = step_dd;
            end

            error_h_log(n2) = error_h;
            error_v_log(n2) = error_h;
            
            heq_00 = heq_00*(1-step*tap_leak) - step*conj(eq_buffer_h)*error_h;
            heq_01 = heq_01*(1-step*tap_leak) - step*conj(eq_buffer_v)*error_h;
            heq_10 = heq_10*(1-step*tap_leak) - step*conj(eq_buffer_h)*error_v;
            heq_11 = heq_11*(1-step*tap_leak) - step*conj(eq_buffer_v)*error_v;
            
            if mod(n2,subsf)==0
                n3 = ceil(n2/subsf);
                heq_00_log(n3,:) = heq_00;
                heq_01_log(n3,:) = heq_01;
                heq_10_log(n3,:) = heq_10;
                heq_11_log(n3,:) = heq_11;
            end
        end
        
    end
    
    if enable_plots_rx
        
        scatterplot(yrx_h_log(end-1e4:end))

        scatterplot(yrx_v_log(end-1e4:end))
        
        figure
        subplot 211
        plot(real(yrx_h_log),'.')
        subplot 212
        plot(real(yrx_v_log),'.')
        
%         figure
%         subplot 211
%         plot(real(error_h_log),'.')
%         subplot 212
%         plot(real(error_h_log),'.')

        figure
        nline = subsf.*(0:length(heq_00_log)-1);
        subplot 221
        plot(nline/BR,abs(heq_00_log))
        subplot 222
        plot(nline/BR,abs(heq_01_log))
        subplot 223
        plot(nline/BR,abs(heq_10_log))
        subplot 224
        plot(nline/BR,abs(heq_11_log))
        
        NFFT = 1024*8;
        fv = 0:fs_rx/NFFT:fs_rx-fs_rx/NFFT; % vector de frecuencias
        HEQ_00 = fft(heq_00,NFFT);
        HEQ_01 = fft(heq_01,NFFT);
        HEQ_10 = fft(heq_10,NFFT);
        HEQ_11 = fft(heq_11,NFFT);
        figure
        subplot 221
        plot(fv/1e9,abs(HEQ_00))
        subplot 222
        plot(fv/1e9,abs(HEQ_01))
        subplot 223
        plot(fv/1e9,abs(HEQ_10))
        subplot 224
        plot(fv/1e9,abs(HEQ_11))
        
    end
    
%     % Filtrado con matched filter
%     rx_h = filter(h_mf, 1, rx_h);
%     rx_v = filter(h_mf, 1, rx_v);
%     
%     % Downsampling
%     rx_h = rx_h(1+down_phase:2:end);
%     rx_v = rx_v(1+down_phase:2:end);
    
    %%
    %---------------%
    % 4D CORRELATOR
    %---------------%

%     rx_h_ref = rx_h;
%     rx_v_ref = rx_v;
    
    if polarization_swap
        ak_hat_h_aux = ak_hat_h;
        ak_hat_v_aux = ak_hat_v;
        ak_hat_h = ak_hat_v_aux;
        ak_hat_v = ak_hat_h_aux;
    end
    
    if polarity_swap_hi
        ak_hat_h = -real(ak_hat_h) + 1j*imag(ak_hat_h);
    end
    
    if polarity_swap_hq
        ak_hat_h = real(ak_hat_h) - 1j*imag(ak_hat_h);
    end
    
    if polarity_swap_vi
        ak_hat_v = -real(ak_hat_v) + 1j*imag(ak_hat_v);
    end
    
    if polarity_swap_vq
        ak_hat_v = real(ak_hat_v) - 1j*imag(ak_hat_v);
    end
    
    if rotation_h ~= 0
        ak_hat_h = ak_hat_h.*exp(1j*rotation_h);
    end
    
    if rotation_v ~= 0
        ak_hat_v = ak_hat_v.*exp(1j*rotation_v);
    end    
    
    LCORR=1000;
    xh_trim = ak_h(end-LCORR+1:end);
    xv_trim = ak_v(end-LCORR+1:end);
    yh_trim = ak_hat_h(end-LCORR+1:end);
    yv_trim = ak_hat_v(end-LCORR+1:end);

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
        ak_hat_h_aux = ak_hat_h;
        ak_hat_v_aux = ak_hat_v;
        ak_hat_h = ak_hat_v_aux;
        ak_hat_v = ak_hat_h_aux;
    end
    
    if M_detection(1,2) ~= 0 || M_detection(1,4) ~= 0 % detect_rotation_h
        ak_hat_h = imag(ak_hat_h) + 1j*real(ak_hat_h);
    end
    
    if M_detection(3,2) ~= 0 || M_detection(3,4) ~= 0 % detect_rotation_v
        ak_hat_v = imag(ak_hat_v) + 1j*real(ak_hat_v);
    end
    
    if M_detection(1,1) < 0 || M_detection(1,2) < 0 || M_detection(1,3) < 0 || M_detection(1,4) < 0 % detect_polarity_swap_hi
        ak_hat_h = -real(ak_hat_h) + 1j*imag(ak_hat_h);
    end
    
    if M_detection(2,1) < 0 || M_detection(2,2) < 0 || M_detection(2,3) < 0 || M_detection(2,4) < 0 % detect_polarity_swap_hq
        ak_hat_h = real(ak_hat_h) - 1j*imag(ak_hat_h);
    end    

    if M_detection(3,1) < 0 || M_detection(3,2) < 0 || M_detection(3,3) < 0 || M_detection(3,4) < 0 % detect_polarity_swap_vi
        ak_hat_v = -real(ak_hat_v) + 1j*imag(ak_hat_v);
    end

    if M_detection(4,1) < 0 || M_detection(4,2) < 0 || M_detection(4,3) < 0 || M_detection(4,4) < 0 % detect_polarity_swap_vq
        ak_hat_v = real(ak_hat_v) - 1j*imag(ak_hat_v);
    end    
    
    if enable_plots_rx
        xh_trim = ak_h(end-LCORR+1:end);
        xv_trim = ak_v(end-LCORR+1:end);
        yh_trim = ak_hat_h(end-LCORR+1:end);
        yv_trim = ak_hat_v(end-LCORR+1:end);

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
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%
    
%     o_data_s.rx_h = rx_h;
%     o_data_s.rx_v = rx_v;
    o_data_s.ak_hat_h = ak_hat_h;
    o_data_s.ak_hat_v = ak_hat_v;
        
end

