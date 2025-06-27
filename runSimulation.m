function runSimulation(rangeStr, ax, tbl, fig, N_users, P_w)
    try
        snr_vals = str2num(rangeStr); %#ok<ST2NM>
        if isempty(snr_vals) || mod(N_users, 2) ~= 0 || P_w <= 0 || P_w >= 1
            error('Invalid inputs');
        end
    catch
        uialert(fig, ['Please ensure:\n' ...
            '- SNR Range is like 0:5:30\n' ...
            '- Number of users is EVEN\n' ...
            '- Power to Weak is between 0 and 1'], 'Input Error');
        return;
    end

    tx_antennas = 4;
    num_bits = 1e4;
    num_iter = 100;
    P_s = 1 - P_w;

    BER_weak_avg = zeros(1, length(snr_vals));
    BER_strong_avg = zeros(1, length(snr_vals));

    for idx = 1:length(snr_vals)
        SNR_dB = snr_vals(idx);
        snr = 10^(SNR_dB / 10);
        total_ber_weak = 0;
        total_ber_strong = 0;

        for run = 1:num_iter
            H = (randn(N_users, tx_antennas) + 1j * randn(N_users, tx_antennas)) / sqrt(2);
            [~, sorted_idx] = sort(vecnorm(H, 2, 2), 'descend');
            clusters = reshape(sorted_idx, 2, []).';

            for c = 1:size(clusters, 1)
                u_s = clusters(c, 1); u_w = clusters(c, 2);
                bits_s = randi([0 1], num_bits, 1); symbols_s = 2 * bits_s - 1;
                bits_w = randi([0 1], num_bits, 1); symbols_w = 2 * bits_w - 1;

                x = [sqrt(P_s)*symbols_s.'; sqrt(P_w)*symbols_w.'];
                H_pair = H([u_s, u_w], :);
                W = (H_pair') / (H_pair * H_pair' + 1e-2 * eye(2));
                W = W ./ vecnorm(W);

                tx = W * x;
                rx = H_pair * tx;
                noise = sqrt(1/(2*snr)) * (randn(size(rx)) + 1j * randn(size(rx)));
                rx_n = rx + noise;

                rx_w = real(rx_n(2, :)); decoded_w = rx_w > 0;
                BER_w = mean(decoded_w ~= bits_w.');

                est_w = sqrt(P_w) * (2 * decoded_w - 1);
                rx_s_clean = real(rx_n(1, :)) - est_w;
                decoded_s = rx_s_clean > 0;
                BER_s = mean(decoded_s ~= bits_s.');

                total_ber_weak = total_ber_weak + BER_w;
                total_ber_strong = total_ber_strong + BER_s;
            end
        end

        BER_weak_avg(idx) = total_ber_weak / (num_iter * size(clusters, 1));
        BER_strong_avg(idx) = total_ber_strong / (num_iter * size(clusters, 1));
        fprintf(' SNR = %2d dB | Weak BER = %.4f | Strong BER = %.4f\n', ...
            SNR_dB, BER_weak_avg(idx), BER_strong_avg(idx));
    end

    % Plot with custom colors and markers
    cla(ax);
    semilogy(ax, snr_vals, BER_weak_avg, '-o', ...
        'LineWidth', 2, 'Color', [0 0.6 1], 'MarkerFaceColor', [0 0.6 1], 'MarkerSize', 8);
    hold(ax, 'on');
    semilogy(ax, snr_vals, BER_strong_avg, '--s', ...
        'LineWidth', 2, 'Color', [1 0.3 0.3], 'MarkerFaceColor', [1 0.3 0.3], 'MarkerSize', 8);
    hold(ax, 'off');

    legend(ax, 'Weak User', 'Strong User', 'Location', 'southwest');
    title(ax, 'BER vs SNR for MIMO-NOMA with SIC');
    xlabel(ax, 'SNR (dB)'); ylabel(ax, 'Bit Error Rate (BER)');
    grid(ax, 'on');

    tbl.Data = [snr_vals(:), BER_weak_avg(:), BER_strong_avg(:)];
end
