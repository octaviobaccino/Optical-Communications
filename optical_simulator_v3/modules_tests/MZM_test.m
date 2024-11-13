clear all
close all


% === Ejemplo 2 ====
% Modelar un MZM IQ
Vpi=6;
swing=0.8; % Se define como el voltaje pico a pico a la entrada de los electrodos del MZI como porcentaje del Vpi. Max swing=2
P0_dBm=20;
Vbias=-Vpi;
ER_inner_dB=inf;
ER_outer_dB=inf;
phase_error=0/180*pi;

% Aux
P0=1e-3*10^(P0_dBm/10);
ER_inner=10^(-ER_inner_dB/20);
ER_outer=10^(-ER_outer_dB/20);

% Genero una señal PAM4 para atacar el MZM
Lsim=100e3;
BR=32e9;
x_rf_i = 2*randi([0 3],Lsim,1)-3;
x_rf_q = 2*randi([0 3],Lsim,1)-3;
x_rf_i =x_rf_i /3*swing/2*Vpi;
x_rf_q =x_rf_q /3*swing/2*Vpi;

P0_inner_I = (1/2+ER_outer)*P0; % Esto es potencia del laser
P0_inner_Q = (1/2-ER_outer)*P0;

% Simulo el MZI de I
p_up = sqrt(1/2+ER_inner)*sqrt(P0_inner_I); % Entrada del modulador de fase UP
p_dn = sqrt(1/2-ER_inner)*sqrt(P0_inner_I); % Entrada del modulador de fase DOWN
y_pm_up = p_up*exp(1j* (x_rf_i+Vbias)*pi/2/Vpi);
y_pm_dn = p_dn*exp(1j*-(x_rf_i+Vbias)*pi/2/Vpi);
y_mzi_i = 1/sqrt(2)*(y_pm_up+y_pm_dn); 

% Simulo el MZI de Q
p_up = sqrt(1/2+ER_inner)*sqrt(P0_inner_Q); % Entrada del modulador de fase UP
p_dn = sqrt(1/2-ER_inner)*sqrt(P0_inner_Q); % Entrada del modulador de fase DOWN
y_pm_up = p_up*exp(1j* (x_rf_q+Vbias)*pi/2/Vpi);
y_pm_dn = p_dn*exp(1j*-(x_rf_q+Vbias)*pi/2/Vpi);
y_mzi_q = 1/sqrt(2)*(y_pm_up+y_pm_dn); 

y_mzi = 1/sqrt(2) * (y_mzi_i + 1j*exp(1j*phase_error)*y_mzi_q); 

% const_theory_pam4=[-0.0929, -0.0929/3, 0.0929/3, 0.0929]; 
const_theory_pam4_i = (2*randi([0 3],1e3,1)-3)/3*0.0929;
const_theory_pam4_q = (2*randi([0 3],1e3,1)-3)/3*0.0929;
figure
plot(const_theory_pam4_i,const_theory_pam4_q,'rx')
hold all
plot(real(y_mzi),imag(y_mzi),'.')

grid on
axis square