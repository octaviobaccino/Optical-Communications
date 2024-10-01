%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
% Created on   : July 2023
% Description  : BER checker
%-----------------------------------------------------------------------------%

function [ber, n_errors] = my_ber_checker(data_v, ref_v, M, guard)
    
    %--------------------------%
    %           ERRORS
    %--------------------------%
    
    if length(data_v) ~= length(ref_v)
        error('The vectors have different length')
    end
    
    %--------------------------%
    %   CONSTANTS & VARIABLES
    %--------------------------%
    
    % Guard
    if guard == 'auto'
        guard = round(length(data_v)/10);
    end
    
    % % CS corrector
    % cfg_cs_s.window_length = 8;
    
    %--------------------------%
    %          PROCESS
    %--------------------------%

    % Align and guard
    alig_delay = finddelay(ref_v, data_v);
    if alig_delay<0
        alig_delay = 0;
    end
    
    ak_guard_v = ref_v(1+guard:end-alig_delay);
    ak_hat_guard_v = data_v(1+guard+alig_delay:end);
    
    % % CS Corrector
    % i_data_s.rx_signal_v = ak_hat_guard_v;
    % i_data_s.tx_signal_v = ak_guard_v;
    % [o_cd_s, ~] = m_cs_corrector(i_data_s, cfg_cs_s);
    % ak_hat_guard_v = o_cd_s.signal_v;
    % ak_guard_v = ak_guard_v(1:length(ak_hat_guard_v));
    
    % QAM to bits
    ak_bit_v = qamdemod(ak_guard_v, M, 'OutputType', 'bit');
    ak_hat_bit_v = qamdemod(ak_hat_guard_v, M, 'OutputType', 'bit');
    
    % BER 
    n_errors = sum(ak_bit_v ~= ak_hat_bit_v);
    ber = n_errors / length(ak_bit_v);
    
end
