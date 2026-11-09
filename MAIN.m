N_mc = 256;                                                 
fc = 7.08e6; c = 3e8; lamda = c / fc;                       
Tr = 0.256;                                                 
Br = 120e3; fs = 1000;                                      
f0 = fc - Br / 2;                                           
Kr = Br / Tr;                                               
t = 0 : 1 / fs : Tr - 1 / fs;                               
N_r = length(t);                                            
Rmax = c / (2 * Br) * N_r;                                  
Vmax = 1 / 4 * lamda / Tr;                                  
Doppler_res = 1 / (Tr * N_mc);                              

N_array = 16; beam_angle = 0;                               
steering_phase = exp(1j * 2*pi *((lamda / 2) / lamda) * (0:N_array - 1)' * sind(beam_angle));

N_target = 20;                                             
R_t = rand(1, N_target) * Rmax;                             
V_t = (rand(1, N_target) * 2 - 1) * Vmax;                   
target_angle = rand(1, N_target) * 5 - 2.5 + beam_angle;    
SNR = -5;                                                   
NoisePower = -SNR + mag2db(sqrt(N_r));
[sr_noise, ~] = func_radar_echo(N_r, Tr, t, N_array, N_target, R_t, V_t, target_angle, zeros(1, N_target), f0, c, lamda, Kr, NoisePower);
sr_noise_dbf = func_LCMV(sr_noise, steering_phase, N_mc, N_r);

[sr_PC] = func_PC(sr_noise_dbf, N_r);
[RD, RD_dB] = func_RD(N_mc, sr_PC);


function [sr_dbf_PC] = func_PC(sr_dbf, N_r)
    range_win = hamming(N_r);
    sr_dbf_windowed = sr_dbf .* range_win';
    sr_dbf_PC = fft(sr_dbf_windowed, N_r, 2);
    sr_dbf_PC = sr_dbf_PC(:, [2:end, 1]);
end


function [mtd, mtd_dB] = func_RD(N_mc, sr_dbf_PC)
    doppler_win = hamming(N_mc);
    sr_dbf_windowed = sr_dbf_PC .* doppler_win;
    mtd = fftshift(fft(sr_dbf_windowed, N_mc, 1), 1);
    mtd = mtd([2:end, 1], :);
    mtd_dB = db(abs(mtd));
end


function sr_noise_dbf = func_LCMV(sr_noise, steering_phase, N_mc, N_r)
    sr_noise_dbf = zeros(N_mc, N_r);
    for kk1 = 1:N_mc
        sr_noise_dbf(kk1, :) = steering_phase' * squeeze(sr_noise(:, kk1, :));
    end
end
