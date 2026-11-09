function [sr_noise, V_t] = func_radar_echo(N_r, Tr, t, N_array, N_target, R_t, V_t1, ...
    target_angle, target_Power, f0, c, lamda, Kr, NoisePower)
    AT = 1; max_doppler = 5;
    
    V_t = zeros(N_target, N_r);
    [b, a] = butter(2, 0.01, 'low');
    for i = 1:N_target
        w = randn(1, N_r);
        x = filter(b, a, w);
        x = (x - mean(x)) * max_doppler * rand();
        V_t(i,:) = V_t1(i) + x;
    end
    
    V_t_shifted = [zeros(N_target,1), V_t(:,1:end-1)];
    V_t_avg = (V_t + V_t_shifted) * 0.5 * Tr;
    R_accumulated = cumsum(V_t_avg, 2);
    
    sr_noise = zeros(N_array, N_r, length(t));
    
    St = AT*exp((1i*2*pi)*(f0*t+Kr/2*t.^2));
    St_conj = conj(St);
    two_pi = 2*pi;
    two_over_c = 2/c;
    phase_base = 1j * two_pi * (lamda/2)/lamda;
    array_indices = (0:N_array-1)';
    target_amp = AT * 10.^(target_Power/20);
    
    phase_shifts = exp(phase_base * array_indices * sind(target_angle));
    
    noise_amp = 10^(NoisePower/20);
    noise_size = [N_array, length(t)];
    
    for i=1:N_r
        sri = zeros(N_array, length(t));
        
        for k=1:N_target
            R_current = R_t(k) + R_accumulated(k, i);
            V_current = V_t(k, i);
            
            tao = two_over_c * (R_current + V_current*t);
            t_delayed = t - tao;
            
            srj = target_amp(k) * exp((1i*two_pi)*(f0*t_delayed + Kr/2*t_delayed.^2));
            
            sri = sri + phase_shifts(:,k) * srj;
        end
        
        noise = normrnd(0, noise_amp, noise_size) + 1j*normrnd(0, noise_amp, noise_size);
        sr_noise(:,i,:) = St_conj .* (sri + noise);
    end
end
