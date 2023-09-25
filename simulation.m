% x = exp(1j*pi*rand(100,1))
% angle(mean(x))
% angle(mean(real(x)) + 1j*mean(imag(x)))
rng(0);
offset = 0.6; % [-pi:0.2:pi];
N = 100000;  % number of measurements
M = 32;     % antenna size
path_ang = -45+90*rand(N,1);
% ensure no wrapping around 2*pi
rel_phase_adj_ant = exp(1j*pi*sin(pi*path_ang/180)); 
% accumulation works (averaging over real and imaginary parts individually)
est_err_adj_ant = angle(mean(real(rel_phase_adj_ant)) + 1j*mean(imag(rel_phase_adj_ant)));
max_err = est_err_adj_ant*(M-1);

Nmeas = [1e3:1e3:1e4];
max_errr_plot = zeros(1,length(Nmeas));
for ii=1:length(Nmeas)
    x = rel_phase_adj_ant(1:Nmeas(ii));
    x2 = angle(mean(real(x)) + 1j*mean(imag(x)));
    max_errr_plot(ii) = x2*(2*M-1);
end
figure; plot(Nmeas, abs(max_errr_plot)); 
xlabel("Number of measurements");
ylabel("Maximum error accumulation (rad)");
title("32x32 URA, 1024 elements");