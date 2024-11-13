function output = laser_generator(p0,LW,delta_f,fs,N,Lsymbs)

% Esta funcion crea un modelo de laser con una potencia promedio
% especifica, offset de portadora y ruido de fase

% p0 = potencia (en watts)
% LW = line width (hasta 500kHz)
% delta_f = offset de portadora (+-2.5GHz)
% fs = frecuencia de muestreo
% Lsymbs = longitud de la simulacion
% N = tasa de sobremuestreo

n = 0:1:(N*Lsymbs-1); % numero total de muestras del laser

p_co = exp(1j*2*pi*delta_f/fs*n);

lw_sigma = sqrt(2*pi*LW*1/fs);
lw_v = cumsum(lw_sigma*randn(1,length(n)));

p_lw = exp(1j*lw_v);

output = sqrt(p0).*p_co.*p_lw;

end

