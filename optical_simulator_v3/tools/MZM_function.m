function [output] = MZM_function(signal,config)
% MZM_FUNCTION Ingreso con una senial compleja, potencia del laser, 
% error de fase, Extintion Ratio y salgo con una senial compleja luego de
% pasar por un modulador de Mach-Zender 

Vpi = config.Vpi;
swing = config.swing; % Se define como el voltaje pico a pico a la entrada de los electrodos del MZI como porcentaje del Vpi. Max swing=2
p0_dbm = config.p0_dbm;
Vbias = -Vpi;
ER_inner_db = config.ER_inner_db;
ER_outer_db = config.ER_outer_db;
phase_error = config.phase_error;

ER_inner=10^(-ER_inner_db/20);
ER_outer=10^(-ER_outer_db/20);

% Laser del MZM

p0 = 10^(p0_dbm/10)*1e-3;

x_rf_i = real(signal);
x_rf_q = imag(signal);
x_rf_i = x_rf_i/max(abs(x_rf_i))*swing*(Vpi/2);
x_rf_q = x_rf_q/max(abs(x_rf_q))*swing*(Vpi/2);

p0_inner_I = (1/2+ER_outer)*p0; % Esto es potencia del laser
p0_inner_Q = (1/2-ER_outer)*p0;

% MZI de I
p_up = sqrt(1/2+ER_inner)*sqrt(p0_inner_I); % Entrada del modulador de fase UP
p_dn = sqrt(1/2-ER_inner)*sqrt(p0_inner_I); % Entrada del modulador de fase DOWN
y_pm_up = p_up*exp(1j* (x_rf_i+Vbias)*pi/(2*Vpi));
y_pm_dn = p_dn*exp(1j*-(x_rf_i+Vbias)*pi/(2*Vpi));
y_mzi_i = 1/sqrt(2)*(y_pm_up+y_pm_dn); 

% MZI de Q
p_up = sqrt(1/2+ER_inner)*sqrt(p0_inner_Q); % Entrada del modulador de fase UP
p_dn = sqrt(1/2-ER_inner)*sqrt(p0_inner_Q); % Entrada del modulador de fase DOWN
y_pm_up = p_up*exp(1j* (x_rf_q+Vbias)*pi/(2*Vpi));
y_pm_dn = p_dn*exp(1j*-(x_rf_q+Vbias)*pi/(2*Vpi));
y_mzi_q = 1/sqrt(2)*(y_pm_up+y_pm_dn); 

% Salida del MZI
output = 1/sqrt(2) * (y_mzi_i + 1j*exp(1j*phase_error)*y_mzi_q); 
end

