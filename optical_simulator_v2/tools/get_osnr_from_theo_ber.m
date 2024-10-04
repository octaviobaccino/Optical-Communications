% Get the OSNR to obtain the theoretical BER

function [osnr_db_v] = get_osnr_from_theo_ber(ber_v, M, BR)

    if M>2
        mod = 'QAM';
    else
        mod = 'PAM';
    end
    
    osnr_db_v = zeros(size(ber_v));
    
    for idx = 1:length(ber_v)
    
        ber = ber_v(idx);
        done = 0;
        ebno_int_db_v = 0:0.1:15;

        % Iteration to obtain ber vector

        while done == 0

           ber_int_v = berawgn(ebno_int_db_v, mod, M);

           if (min(ber_int_v) < ber) && (max(ber_int_v) > ber)

               done = 1;

           elseif min(ber_int_v) > ber

               ebno_int_db_v = ebno_int_db_v + 0.1;

           elseif max(ber_int_v) < ber    

               ebno_int_db_v = ebno_int_db_v - 0.1;

           end

        end

        % Interpolation
        ebno_db = interp1(log10(ber_int_v), ebno_int_db_v, log10(ber));
        snr_db = ebno_db + 10*log10(log2(M));
        osnr_db = snr_db + 10*log10(BR/12.5e9);
        osnr_db_v(idx) = osnr_db;

    end
end
