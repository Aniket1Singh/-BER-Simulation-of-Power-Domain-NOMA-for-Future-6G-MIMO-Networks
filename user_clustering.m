function clusters = user_clustering(H)
    % Clusters users based on channel strength (norm of rows)
    [~, idx] = sort(vecnorm(H, 2, 2));  % Sort users by channel strength
    clusters = reshape(idx, [], 2);     % Pair weak+strong
end
