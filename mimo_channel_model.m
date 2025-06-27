function H = mimo_channel_model(N_users, tx_antennas)
    % Generates Rayleigh fading MIMO channel matrix
    H = (randn(N_users, tx_antennas) + 1j * randn(N_users, tx_antennas)) / sqrt(2);
end
