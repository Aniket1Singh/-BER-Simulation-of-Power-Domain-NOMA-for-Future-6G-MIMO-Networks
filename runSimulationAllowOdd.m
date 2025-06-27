function runSimulationAllowOdd(snrStr, ax, tbl, users, pw)
    try
        snr_vals = eval(snrStr);
        if isempty(snr_vals)
            error('SNR range is empty.');
        end
    catch
        uialert(gcf, 'Invalid SNR range format. Use format like 0:5:30', 'Input Error');
        return;
    end

    % Basic parameters
    nIter = 1000;
    nSymbols = 1e4;
    weakBER = zeros(size(snr_vals));
    strongBER = zeros(size(snr_vals));

    for i = 1:length(snr_vals)
        snr = snr_vals(i);
        noiseVar = 10^(-snr/10);

        weakErrors = 0;
        strongErrors = 0;

        for k = 1:nIter
            % Simulate random bits for all users
            data = randi([0 1], users, nSymbols);

            % Superposition coding
            weakData = data(1,:);
            strongData = data(2,:);
            txSig = sqrt(pw)*bpsk(weakData) + sqrt(1-pw)*bpsk(strongData);

            % AWGN channel
            rx = txSig + sqrt(noiseVar/2)*randn(1,nSymbols);

            % SIC at strong user
            estWeak = bpsk(rx > 0);  % weak decoded
            rxCleaned = rx - sqrt(pw)*bpsk(estWeak);
            estStrong = bpsk(rxCleaned > 0);

            % BER
            weakErrors = weakErrors + sum(estWeak ~= weakData);
            strongErrors = strongErrors + sum(estStrong ~= strongData);
        end

        weakBER(i) = weakErrors / (nSymbols * nIter);
        strongBER(i) = strongErrors / (nSymbols * nIter);
    end

    % Plot
    cla(ax);
    semilogy(ax, snr_vals, weakBER, 'ro-', ...
             snr_vals, strongBER, 'bs--', 'LineWidth', 1.5);
    legend(ax, 'Weak User', 'Strong User');
    title(ax, 'BER vs SNR for MIMO-NOMA');
    xlabel(ax, 'SNR (dB)');
    ylabel(ax, 'BER');
    grid(ax, 'on');

    % Update table
    tbl.Data = num2cell([snr_vals(:), weakBER(:), strongBER(:)]);
end

function out = bpsk(bits)
    out = 2*bits - 1;
end
