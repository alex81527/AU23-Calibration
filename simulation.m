rng(4);
offset = 0.6; % [-pi:0.2:pi];
N = 100000;   % number of measurements
URA_size = [36, 64, 256, 1024, 4096];
path_ang = -45+90*rand(N,1);
% ensure no wrapping around 2*pi (slides 28)
rel_phase_adj_ant = exp(1j*pi*sin(pi*path_ang/180)); 
% Per antenna, accumulating measurements into one complex number works 
est_err_adj_ant = angle(mean(rel_phase_adj_ant));

M = sqrt(URA_size);     % antenna size, per side
Nmeas = [1e2:1e2:1e4];
max_accumulation_err = zeros(length(Nmeas), length(URA_size));
% max_accumulation_err = est_err_adj_ant*(2*(M-1));

for ii=1:length(Nmeas)
    x = angle(mean(rel_phase_adj_ant(1:Nmeas(ii))));
    max_accumulation_err(ii, :) = x*(2*(ceil(M/2)-1));
end

figure; plot(Nmeas, abs(max_accumulation_err)); 
xlabel("Number of measurements");
ylabel("Maximum error accumulation (rad)");
legend("URA size "+string(URA_size));
title("One realization of random measurements");