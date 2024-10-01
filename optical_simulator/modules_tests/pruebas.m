close all;

config.BR = 64e9;
BR = config.BR;
config.config_rx.enable_plots_rx = 0;
config.config_ch.osnr_db = 18;

% tx_odata = optical_TX(config);
% 
% i_channel_s.tx_out_h = tx_odata.tx_out_h;
% i_channel_s.tx_out_v = tx_odata.tx_out_v;
% 
% ch_odata = optical_CH(i_channel_s,config);
% 
% rx_odata = optical_RX(ch_odata,config);
% 
% [pxx1,f1] = pwelch(i_channel_s.tx_out_h,hanning(1024*8/2),0,1024*8,2*BR);
% [pxx2,f2] = pwelch(ch_odata.xup_ch_h,hanning(1024*8/2*2),0,1024*8,2*2*BR);
% [pxx3,f3] = pwelch(ch_odata.yup_ch_h,hanning(1024*8/2*2),0,1024*8,2*2*BR);
% [pxx4,f4] = pwelch(ch_odata.ch_out_h,hanning(1024*8/2*2),0,1024*8,2*2*BR);
% 
% 
% figure
% plot(f1,pxx1)
% hold all
% plot(f3,pxx3,'--')
% 
% figure
% plot(f1,pxx1)
% hold all
% plot(f4,pxx4,'--')

% figure
% plot(f1,pxx1)
% hold all
% plot(f2,pxx2,'--')

odata = optical_simulator_v1(config);

