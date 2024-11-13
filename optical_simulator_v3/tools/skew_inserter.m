function [odata] = skew_inserter(idata, skew, fs)

% idata es una senial compleja con contribucion de I y Q
% skew esta en segundos
% fs es la fcia de muestreo del canal

tx_ui = skew*fs;
hsk_i = my_rcosine(1, 1, 0.1, 100, tx_ui/2);
hsk_q = my_rcosine(1, 1, 0.1, 100, -tx_ui/2);

% figure
% timeline = (-length(hsk_i)/2:length(hsk_i)/2-1)./fs;
% plot(timeline/1e-12, hsk_i); hold all
% plot(timeline/1e-12, hsk_q); 

skew_gen_input_i = real(idata);
skew_gen_input_q = imag(idata);
skew_gen_output_i = filter(hsk_i, 1, skew_gen_input_i);
skew_gen_output_q = filter(hsk_q, 1, skew_gen_input_q);
odata = skew_gen_output_i + 1j*skew_gen_output_q;
end