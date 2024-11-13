%% Gain Imbalance sweep processing

clear; close all;

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

file = [out_dir, 'cfg.mat'];
load(file);

n_ER = length(ER_outer_db_v);
n_ber = length(theo_ber_v);

M = config.M;
BR = config.BR;

ber_est_m = zeros(n_ER,n_ber);
ber_theo_m = zeros(n_ER,n_ber);

%% Read data

for idx_ER = 1:n_ER
    
    ER_outer_db = ER_outer_db_v(idx_ER);
    
    name = sprintf('out_ER_outer_db_%2.2f', ER_outer_db);
    file = [out_dir, name, '.mat'];
    load(file);
            
    for idx_ber = 1:n_ber

        ber_est_m(idx_ER,idx_ber) = out_c{idx_ber}.ber_est;
        ber_theo_m(idx_ER,idx_ber) = out_c{idx_ber}.ber_theo;

    end
    
end

%% Plots

% SNR Penalty vs ER_outer_db

ber_int = 1e-2;

snr_loss_db_v = zeros(n_ER,1);

ber_theo_v = zeros(n_ber,1);
ber_est_v = zeros(n_ber,1);


% Procesamiento

    for idx_ER = 1:n_ER

        ber_theo_v = ber_theo_m(idx_ER,:);
        ber_est_v = ber_est_m(idx_ER,:);
        
        osnr_db_v = get_osnr_from_theo_ber(ber_theo_v,M,BR);

        osnr_sim_db = interp1(log10(ber_est_v), osnr_db_v, log10(ber_int));
        osnr_theo_db = interp1(log10(ber_theo_v), osnr_db_v, log10(ber_int));

        snr_loss_db_v(idx_ER) =  osnr_sim_db - osnr_theo_db;

    end

fz = 15;
idx_leg = 1;
color_c = { [0 0.45 0.74]
            [0.85 0.33 0.1]
            [0.93 0.69 0.13]
            [0.49 0.18 0.56]
            [0.47 0.67 0.19]
            [0.30 0.75 0.93]
            [0.64 0.08 0.18]};
        
figure
p = plot(ER_outer_db_v, snr_loss_db_v, '-o', 'Linewidth', 1.1);
p.MarkerFaceColor = color_c{idx_leg};
p.MarkerEdgeColor = 'k';
p.Color = color_c{idx_leg};
hold on; grid on;
xlabel('ER\_outer [dB]', 'Interpreter','latex','FontSize', fz);
ylab = sprintf('ROSNR [dB] @ BER = %.1e',ber_int);
ylabel(ylab, 'Interpreter','latex','FontSize', fz);
tit = ['OSNR loss vs Extintion Ratio. ',sprintf('BR=%d[GBd]', BR/1e9)];
title(tit,'Interpreter','latex','FontSize', fz);
out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
save_file = [out_dir 'gain_imbalance.png'];
saveas(gcf,save_file);
    