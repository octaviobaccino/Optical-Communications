%% Step sweep processing

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

n_BR = length(BR_v);
n_steps = length(step_dd_v);
n_ber = length(theo_ber_v);

M = config.M;

BR_cell = cell(n_BR,1);

ber_est_m = zeros(n_steps,n_ber);
ber_theo_m = zeros(n_steps,n_ber);

%% Read data

for idx_BR = 1:n_BR
    
    BR = BR_v(idx_BR);
    
    out_dir_BR = [out_dir(1:end-1), '\'];
    name = sprintf('out_BR_%d/', BR/1e9);
    out_dir_BR = [out_dir_BR, name];
    
    if ~exist(out_dir_BR,'dir')
        mkdir(out_dir_BR);
    end
    
    for idx_step = 1:n_steps
        
        step_dd = step_dd_v(idx_step);
        
        name = sprintf('out_step_dd_%e',step_dd);
        file = [out_dir_BR, name, '.mat'];
        load(file);
        
        for idx_ber = 1:n_ber
            
            ber_est_m(idx_step,idx_ber) = out_c{idx_ber}.ber_est;
            ber_theo_m(idx_step,idx_ber) = out_c{idx_ber}.ber_theo;
            
        end
        
        aux_s.ber_est_m = ber_est_m;
        aux_s.ber_theo_m = ber_theo_m;
        
    end
    
    BR_cell{idx_BR,1} = aux_s; 
    
end

%% Plots

fz = 15;
idx_leg = 1;
color_c = { [0 0.45 0.74]
            [0.85 0.33 0.1]
            [0.93 0.69 0.13]
            [0.49 0.18 0.56]
            [0.47 0.67 0.19]
            [0.30 0.75 0.93]
            [0.64 0.08 0.18]};

for idx_BR = 1:n_BR
    
    BR = BR_v(idx_BR);
    
    figure

    osnr_db_v = get_osnr_from_theo_ber(theo_ber_v,M,BR);

    ber_theo_v = BR_cell{idx_BR,1}.ber_theo_m(idx_step,:);
    
    semilogy(osnr_db_v, ber_theo_v)
    hold all;
    grid on;
    
    for idx_step = 1:n_steps
        ber_est_v = BR_cell{idx_BR,1}.ber_est_m(idx_step,:);
        semilogy(osnr_db_v, ber_est_v)
    end
end
    
