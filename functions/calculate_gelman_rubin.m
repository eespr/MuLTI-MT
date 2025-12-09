function [R_hat] = calculate_gelman_rubin(chain_data, burn_in)
% [R_hat] = calculate_gelman_rubin(chain_data, burn_in)
%
% Calculates the Gelman-Rubin statistic (R_hat), also known as the
% potential scale reduction factor (PSRF), to assess MCMC convergence.
%
% INPUTS:
%   chain_data = (M x N) matrix, where M is the number of chains
%                and N is the total number of samples per chain.
%   burn_in    = Scalar, the number of samples to discard as burn-in.
%
% OUTPUT:
%   R_hat      = The Gelman-Rubin statistic. A value < 1.1 (or 1.05)
%                is typically considered to indicate convergence.
%
% Based on the method described by Gelman and Rubin (1992) and
% Brooks and Gelman (1998).

% Select post-burn-in data
if burn_in >= size(chain_data, 2)
    error('Burn-in is greater than or equal to the number of samples.');
end
post_burn_in_data = chain_data(:, burn_in+1:end);

[M, N] = size(post_burn_in_data); % M = number of chains, N = number of post-burn-in samples

if M <= 1
    warning('Gelman-Rubin statistic requires at least 2 chains. Returning R_hat = NaN.');
    R_hat = NaN;
    return;
end

% 1. Calculate the mean of each chain
chain_means = mean(post_burn_in_data, 2, 'omitnan');

% 2. Calculate the grand mean (mean of all samples from all chains)
grand_mean = mean(chain_means, 'omitnan');

% 3. Calculate the Between-chain variance (B)
B = (N / (M - 1)) * sum((chain_means - grand_mean).^2, 'omitnan');

% 4. Calculate the Within-chain variance (W)
% First, get the variance for each chain
chain_variances = var(post_burn_in_data, 0, 2, 'omitnan'); % 0 = normalize by N-1
% Then, average them
W = mean(chain_variances, 'omitnan');

% Handle cases where variance is zero (e.g., all chains stuck at one value)
if W == 0
    % If W is 0, chains haven't explored.
    % If B is also 0, they are all stuck at the same point (R_hat=1 is ok)
    % If B > 0, they are stuck at different points (R_hat=Inf is appropriate)
    if B == 0
        R_hat = 1;
    else
        R_hat = Inf;
    end
    return;
end

% 5. Estimate the pooled posterior variance (V_hat)
V_hat = ((N - 1) / N) * W + (B / N);

% 6. Calculate the Potential Scale Reduction Factor (PSRF or R_hat)
R_hat = sqrt(V_hat / W);

% Handle potential numerical issues
if isnan(R_hat)
    R_hat = 1;
end

end