clc; clear; close all;

%% System Parameters
N_users = 8;               % Total users
tx_antennas = 4;           % Number of transmit antennas (BS)
rx_antennas = 1;           % 1 antenna per user
num_bits = 1e4;            % Bits per user
mod_order = 2;             % BPSK
SNR_dB_range = 0:5:30;     % SNR sweep range in dB
num_runs = 100;            % Monte Carlo iterations per SNR

%% Initialize Result Storage
BER_weak_avg = zeros(size(SNR_dB_range));
BER_strong_avg = zeros(size(SNR_dB_range));

%% SNR Sweep
for snr_idx = 1:length(SNR_dB_range)
    SNR_dB = SNR_dB_range(snr_idx);
    snr = 10^(SNR_dB/10);

    total_BER_weak = 0;
    total_BER_strong = 0;
    total_clusters = 0;

    for run = 1:num_runs
        %% Generate Random Channel Matrix
        H = (randn(N_users, tx_antennas) + 1j * randn(N_users, tx_antennas)) / sqrt(2);

        %% User Clustering (based on channel norm)
        [~, sorted_idx] = sort(vecnorm(H, 2, 2), 'descend');
        clusters = reshape(sorted_idx, 2, []).';

        %% Loop Through Clusters
        for c = 1:size(clusters, 1)
            u_strong = clusters(c, 1);
            u_weak   = clusters(c, 2);

            % Generate Data
            bits_s = randi([0 1], num_bits, 1); symbols_s = 2 * bits_s - 1;
            bits_w = randi([0 1], num_bits, 1); symbols_w = 2 * bits_w - 1;

            % Power allocation
            P_s = 0.3; P_w = 0.7;
            x_combined = [sqrt(P_s)*symbols_s.'; sqrt(P_w)*symbols_w.'];  % 2 x num_bits

            % Channel matrix for cluster
            H_pair = H([u_strong, u_weak], :);

            % Regularized Zero-Forcing Beamforming
            W = (H_pair') / (H_pair * H_pair' + 1e-2 * eye(2));
            W = W ./ vecnorm(W);  % Normalize columns

            % Transmit signal
            tx_bf = W * x_combined;     % [tx_ant x num_bits]
            rx_cluster = H_pair * tx_bf; % [2 x num_bits]

            % Add AWGN
            noise = sqrt(1/(2*snr)) * (randn(size(rx_cluster)) + 1j * randn(size(rx_cluster)));
            rx_noisy = rx_cluster + noise;

            % Weak user decoding
            rx_w = real(rx_noisy(2, :));
            decoded_w = rx_w > 0;
            BER_w = mean(decoded_w ~= bits_w.');

            % Strong user decoding (SIC)
            est_w = sqrt(P_w) * (2 * decoded_w - 1);
            rx_s_clean = real(rx_noisy(1, :)) - est_w;
            decoded_s = rx_s_clean > 0;
            BER_s = mean(decoded_s ~= bits_s.');

            % Accumulate results
            total_BER_weak = total_BER_weak + BER_w;
            total_BER_strong = total_BER_strong + BER_s;
            total_clusters = total_clusters + 1;
        end
    end

    % Average BER over all clusters and runs
    BER_weak_avg(snr_idx) = total_BER_weak / total_clusters;
    BER_strong_avg(snr_idx) = total_BER_strong / total_clusters;

    fprintf("Done for SNR = %d dB: Weak BER = %.4f, Strong BER = %.4f\n", ...
        SNR_dB, BER_weak_avg(snr_idx), BER_strong_avg(snr_idx));
end

%% Plotting
figure;
semilogy(SNR_dB_range, BER_weak_avg, '-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB_range, BER_strong_avg, '-x', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for MIMO-NOMA with SIC');
legend('Weak User', 'Strong User');
grid on;
