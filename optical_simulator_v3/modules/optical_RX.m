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
    config.step_cma = 2^-9;
    config.step_dd = 2^-6;
    
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
    
%     h_mf = root_raised_cosine(BR/2, fs_rx, rolloff, pulse_shaping_ntaps, 0);
    
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
        
        if n < fix(Lrx/5)
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
        grid on
        title('Constelacion a la entrada del demapper','Interpreter','Latex')
        
        scatterplot(yrx_v_log(end-1e4:end));
        grid on
        title('Constelacion a la entrada del demapper','Interpreter','Latex')
        
        
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
        plot(nline/BR/1e-6,abs(heq_00_log))
        grid on;
%         xlabel('Time[us]', 'Interpreter','latex','FontSize', 15);
        ylabel('Taps amplitude', 'Interpreter','latex','FontSize', 15);
        subplot 222
        plot(nline/BR/1e-6,abs(heq_01_log))
        grid on;
%         xlabel('Time[us]', 'Interpreter','latex','FontSize', 15);
        subplot 223
        plot(nline/BR/1e-6,abs(heq_10_log))
        grid on;
        xlabel('Time[us]', 'Interpreter','latex','FontSize', 15);
        ylabel('Taps amplitude', 'Interpreter','latex','FontSize', 15);
        subplot 224
        plot(nline/BR/1e-6,abs(heq_11_log))
        grid on;
        xlabel('Time[us]', 'Interpreter','latex','FontSize', 15);
        
        tit = ['Equlizaer taps. ',sprintf('BR=%d[GBd]', BR/1e9) '. fsop=200[KHz]'];
        sgtitle(tit,'Interpreter','latex','FontSize', 15);

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
    
    % Correlador 4 dimensiones
    
    LCORR=1000;
    [ak_hat_h,ak_hat_v] = correlator_4D(ak_h,ak_v,ak_hat_h,ak_hat_v,LCORR,enable_plots_rx);
    
    %--------------------------%
    %         OUTPUT
    %--------------------------%
    
%     o_data_s.rx_h = rx_h;
%     o_data_s.rx_v = rx_v;
    o_data_s.ak_hat_h = ak_hat_h;
    o_data_s.ak_hat_v = ak_hat_v;
        
end

