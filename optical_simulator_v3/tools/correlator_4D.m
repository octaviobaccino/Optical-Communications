function [output1,output2] = correlator_4D(ak_h,ak_v,ak_hat_h,ak_hat_v,LCORR,enable_plots)
% 4D_CORRELATOR enter with signals of 2 polarizations (h and v)
% transmitted and received symbols and the lenght of polarization

% Tomo solamente la ultima parte de la senial
xh_trim = ak_h(end-LCORR+1:end);
xv_trim = ak_v(end-LCORR+1:end);
yh_trim = ak_hat_h(end-LCORR+1:end);
yv_trim = ak_hat_v(end-LCORR+1:end);

% Hago las correlaciones de todas las posibles mezclas de senial
corr_hi_hi = xcorr(real(xh_trim), real(yh_trim))./var(real(xh_trim))./LCORR;
corr_hi_hq = xcorr(real(xh_trim), imag(yh_trim))./var(real(xh_trim))./LCORR;
corr_hi_vi = xcorr(real(xh_trim), real(yv_trim))./var(real(xh_trim))./LCORR;
corr_hi_vq = xcorr(real(xh_trim), imag(yv_trim))./var(real(xh_trim))./LCORR;

corr_hq_hi = xcorr(imag(xh_trim), real(yh_trim))./var(imag(xh_trim))./LCORR;
corr_hq_hq = xcorr(imag(xh_trim), imag(yh_trim))./var(imag(xh_trim))./LCORR;
corr_hq_vi = xcorr(imag(xh_trim), real(yv_trim))./var(imag(xh_trim))./LCORR;
corr_hq_vq = xcorr(imag(xh_trim), imag(yv_trim))./var(imag(xh_trim))./LCORR;

corr_vi_hi = xcorr(real(xv_trim), real(yh_trim))./var(real(xv_trim))./LCORR;
corr_vi_hq = xcorr(real(xv_trim), imag(yh_trim))./var(real(xv_trim))./LCORR;
corr_vi_vi = xcorr(real(xv_trim), real(yv_trim))./var(real(xv_trim))./LCORR;
corr_vi_vq = xcorr(real(xv_trim), imag(yv_trim))./var(real(xv_trim))./LCORR;

corr_vq_hi = xcorr(imag(xv_trim), real(yh_trim))./var(imag(xv_trim))./LCORR;
corr_vq_hq = xcorr(imag(xv_trim), imag(yh_trim))./var(imag(xv_trim))./LCORR;
corr_vq_vi = xcorr(imag(xv_trim), real(yv_trim))./var(imag(xv_trim))./LCORR;
corr_vq_vq = xcorr(imag(xv_trim), imag(yv_trim))./var(imag(xv_trim))./LCORR;

% Plots de las correlaciones sin arreglar
if enable_plots
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

% Tomo los maximos de cada correlacion y su posicion
[max_hi_hi,pos_hi_hi] = max(abs(corr_hi_hi));
[max_hi_hq,pos_hi_hq] = max(abs(corr_hi_hq));
[max_hi_vi,pos_hi_vi] = max(abs(corr_hi_vi));
[max_hi_vq,pos_hi_vq] = max(abs(corr_hi_vq));

[max_hq_hi,pos_hq_hi] = max(abs(corr_hq_hi));
[max_hq_hq,pos_hq_hq] = max(abs(corr_hq_hq));
[max_hq_vi,pos_hq_vi] = max(abs(corr_hq_vi));
[max_hq_vq,pos_hq_vq] = max(abs(corr_hq_vq));

[max_vi_hi,pos_vi_hi] = max(abs(corr_vi_hi));
[max_vi_hq,pos_vi_hq] = max(abs(corr_vi_hq));
[max_vi_vi,pos_vi_vi] = max(abs(corr_vi_vi));
[max_vi_vq,pos_vi_vq] = max(abs(corr_vi_vq));

[max_vq_hi,pos_vq_hi] = max(abs(corr_vq_hi));
[max_vq_hq,pos_vq_hq] = max(abs(corr_vq_hq));
[max_vq_vi,pos_vq_vi] = max(abs(corr_vq_vi));
[max_vq_vq,pos_vq_vq] = max(abs(corr_vq_vq));

% Creo una matriz para con los maximos y su signo
M_corr = zeros(4,4);

M_corr(1,1) = (max_hi_hi>0.5)*sign(corr_hi_hi(pos_hi_hi));
M_corr(1,2) = (max_hi_hq>0.5)*sign(corr_hi_hq(pos_hi_hq));
M_corr(1,3) = (max_hi_vi>0.5)*sign(corr_hi_vi(pos_hi_vi));
M_corr(1,4) = (max_hi_vq>0.5)*sign(corr_hi_vq(pos_hi_vq));

M_corr(2,1) = (max_hq_hi>0.5)*sign(corr_hq_hi(pos_hq_hi));
M_corr(2,2) = (max_hq_hq>0.5)*sign(corr_hq_hq(pos_hq_hq));
M_corr(2,3) = (max_hq_vi>0.5)*sign(corr_hq_vi(pos_hq_vi));
M_corr(2,4) = (max_hq_vq>0.5)*sign(corr_hq_vq(pos_hq_vq));

M_corr(3,1) = (max_vi_hi>0.5)*sign(corr_vi_hi(pos_vi_hi));
M_corr(3,2) = (max_vi_hq>0.5)*sign(corr_vi_hq(pos_vi_hq));
M_corr(3,3) = (max_vi_vi>0.5)*sign(corr_vi_vi(pos_vi_vi));
M_corr(3,4) = (max_vi_vq>0.5)*sign(corr_vi_vq(pos_vi_vq));

M_corr(4,1) = (max_vq_hi>0.5)*sign(corr_vq_hi(pos_vq_hi));
M_corr(4,2) = (max_vq_hq>0.5)*sign(corr_vq_hq(pos_vq_hq));
M_corr(4,3) = (max_vq_vi>0.5)*sign(corr_vq_vi(pos_vq_vi));
M_corr(4,4) = (max_vq_vq>0.5)*sign(corr_vq_vq(pos_vq_vq));

% Multiplicando la matriz por la senial recibida por la matriz recupero la
% la senial original
ak_hat_hi = real(ak_hat_h)*M_corr(1,1) + imag(ak_hat_h)*M_corr(1,2)...
    + real(ak_hat_v)*M_corr(1,3) + imag(ak_hat_v)*M_corr(1,4);

ak_hat_hq = real(ak_hat_h)*M_corr(2,1) + imag(ak_hat_h)*M_corr(2,2)...
    + real(ak_hat_v)*M_corr(2,3) + imag(ak_hat_v)*M_corr(2,4);

ak_hat_vi = real(ak_hat_h)*M_corr(3,1) + imag(ak_hat_h)*M_corr(3,2)...
    + real(ak_hat_v)*M_corr(3,3) + imag(ak_hat_v)*M_corr(3,4);    

ak_hat_vq = real(ak_hat_h)*M_corr(4,1) + imag(ak_hat_h)*M_corr(4,2)...
    + real(ak_hat_v)*M_corr(4,3) + imag(ak_hat_v)*M_corr(4,4);

ak_hat_h = ak_hat_hi + 1j*ak_hat_hq;
ak_hat_v = ak_hat_vi + 1j*ak_hat_vq;

% plots de la senial arreglada
if enable_plots
    xh_trim = ak_h(end-LCORR+1:end);
    xv_trim = ak_v(end-LCORR+1:end);
    yh_trim = ak_hat_h(end-LCORR+1:end);
    yv_trim = ak_hat_v(end-LCORR+1:end);

    corr_hi_hi = xcorr(real(xh_trim), real(yh_trim))./var(real(xh_trim))./LCORR;
    corr_hi_hq = xcorr(real(xh_trim), imag(yh_trim))./var(real(xh_trim))./LCORR;
    corr_hi_vi = xcorr(real(xh_trim), real(yv_trim))./var(real(xh_trim))./LCORR;
    corr_hi_vq = xcorr(real(xh_trim), imag(yv_trim))./var(real(xh_trim))./LCORR;

    corr_hq_hi = xcorr(imag(xh_trim), real(yh_trim))./var(imag(xh_trim))./LCORR;
    corr_hq_hq = xcorr(imag(xh_trim), imag(yh_trim))./var(imag(xh_trim))./LCORR;
    corr_hq_vi = xcorr(imag(xh_trim), real(yv_trim))./var(imag(xh_trim))./LCORR;
    corr_hq_vq = xcorr(imag(xh_trim), imag(yv_trim))./var(imag(xh_trim))./LCORR;

    corr_vi_hi = xcorr(real(xv_trim), real(yh_trim))./var(real(xv_trim))./LCORR;
    corr_vi_hq = xcorr(real(xv_trim), imag(yh_trim))./var(real(xv_trim))./LCORR;
    corr_vi_vi = xcorr(real(xv_trim), real(yv_trim))./var(real(xv_trim))./LCORR;
    corr_vi_vq = xcorr(real(xv_trim), imag(yv_trim))./var(real(xv_trim))./LCORR;

    corr_vq_hi = xcorr(imag(xv_trim), real(yh_trim))./var(imag(xv_trim))./LCORR;
    corr_vq_hq = xcorr(imag(xv_trim), imag(yh_trim))./var(imag(xv_trim))./LCORR;
    corr_vq_vi = xcorr(imag(xv_trim), real(yv_trim))./var(imag(xv_trim))./LCORR;
    corr_vq_vq = xcorr(imag(xv_trim), imag(yv_trim))./var(imag(xv_trim))./LCORR;

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
    
output1 = ak_hat_h;
output2 = ak_hat_v;
    
end

