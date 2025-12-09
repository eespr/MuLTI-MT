%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  MuLTI MT                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab script runs a transdimensional MCMC 
% inversion for 1-D MT data 

% Author: Siobhan Killingbeck, School of Biosciences, Geography and Physics
% Swansea University (2025)

% It is based on Matlab code MuLTI written by:
% Siobhan Killingbeck and Phil Livermore 
% School of Earth and Environemnt, The University of Leeds

% The physical model consists of internal layers, each with an associated 
% resistivity (R), sampled as log(R) defined by Voronoi nuclei.

% The domain is divided into a number of layers, num_layers (which could be one)
% each with its own prior distribution on R.

% Each layer has a special nuclei that cannot leave its layer ("confined" nuclei).
% npt is the number of "floating" nuclei that can change layer.

% Inverts using COMPLEX IMPEDANCE (Z_real, Z_imag) but plots AR/Ph.
% Inverts the determinant (Zssq), ZXY or ZYX 

%==========================================================================
% core function: run_multiple_chains_Zinvert.m 
%==========================================================================

clear all
close all

%% %%%%%%%%%%%%%%%%%%%%%%% ADD PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add data paths and functions paths
% addpath(genpath('functions')
load('cmap_blue.mat') % <<< Using cm_blue white to blue

%% %%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Z_real, Z_imag, and their errors for inversion
[d, units_code] = load_data_edi('INSERT_FILE_NAME_HERE.edi');

ZZ = 1; % ZZ=1=Zqq / ZZ=2=ZXY /ZZ=3=ZYX

error_floor_percentage = 0.05; % 5% error floor (adjust if needed)

% ------------------------------------------------------------------------
frequency = d.f(:); 
Zxx = d.Z(:, 1); Zxy = d.Z(:, 2); Zyx = d.Z(:, 3); Zyy = d.Z(:, 4);
Zerr_xx = d.Zerr(:, 1); Zerr_xy = d.Zerr(:, 2); Zerr_yx = d.Zerr(:, 3); Zerr_yy = d.Zerr(:, 4);

if ZZ == 1 %Zssq
    % CALCULATE DETERMINANT Zssq
    Z_data = (sqrt((Zxx.^2 + Zxy.^2 + Zyx.^2 + Zyy.^2)./2));
    % Calculate the approximate single error value dZ 
    dZ = real(sqrt(0.5*sqrt(Zerr_xx.^4 + Zerr_xy.^4 + Zerr_yx.^4 + Zyy.^4)));
elseif ZZ == 2 % ZXY
    Z_data = Zxy;
    % Calculate the approximate single error value dZ 
    dZ = Zerr_xy;
elseif ZZ == 3 % ZYX
    Z_data = -Zyx;
    % Calculate the approximate single error value dZ 
    dZ = Zerr_yx;
end

% --- Data to keep for FINAL PLOTTING (AR/Ph space) ---
mua = 4*pi*10^-7;
w = 2*pi*frequency;
dataAR_plot = (1./(w*mua)).*abs(Z_data).^2; % Apparent Resistivity for plotting
dataPh_plot = (atan2(imag(Z_data),real(Z_data))*(180/pi)); % Phase for plotting
rhoaerr_plot = (2*dataAR_plot.*dZ)./abs(Z_data); % AR error for plotting
phierr_plot = (180/pi)*(dZ)./abs(Z_data);      % Phase error for plotting

% --- Data FOR INVERSION (Real/Imag Z space) ---
% Assign approximate error dZ to both real and imaginary
dZ_real = dZ; 
dZ_imag = dZ; 

% Apply error floor
Z_magnitude = abs(Z_data);
Z_error_floor = Z_magnitude * error_floor_percentage;
weightingZ_real = max(dZ_real, Z_error_floor);
weightingZ_imag = max(dZ_imag, Z_error_floor);

nd = numel(dataAR_plot); % Still based on number of frequencies

% --- Package data for the INVERSION function ---
data.frequency = frequency;
data.Z_data = Z_data; % Pass complex Z
data.weightingZ_real = weightingZ_real; % Pass real error/weight
data.weightingZ_imag = weightingZ_imag; % Pass imag error/weight
data.nd = nd;
data.running_mode = 1; % 1 -> find posterior; 0 -> Find priors.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     DEFINE PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_chains = 10; % Define number of chains to run

params.burn_in = 10000; % burn-in period
params.nsample = 1000000; % total number of samples

% Uniform prior on depths for nuclei
priors.depth_min=0; 
priors.depth_max=100000; % in meters
params.dis = 500; % steps to discretize the model.
priors.npt_max = 100; % max number of floating nuclei
priors.npt_min = 0; % min number of floating nuclei

% Depth constraining parameter
% Unconstrained: priors.num_layers = 1 
% Constrained: priors.num_layers > 1
priors.num_layers = 1;  

priors.layer_depths = zeros(priors.num_layers-1,1);
priors.Rmin = zeros(priors.num_layers,1);
priors.Rmax = zeros(priors.num_layers,1);

% Define the priors and layer geometry
if priors.num_layers == 3
    priors.layer_depths(1) = 3388 ; % layer 1 depth.
    priors.layer_depths(2) = NaN ; % layer 2 depth.
    priors.Rmin(1) = log10(5000); priors.Rmax(1) = log10(1000000); 
    priors.Rmin(2) = log10(0.1); priors.Rmax(2) = log10(1000);
    priors.Rmin(3) = log10(0.1); priors.Rmax(3) = log10(1000000);   
elseif priors.num_layers == 2
    priors.layer_depths(1) = 3388 ; % layer 1 depth.
    priors.Rmin(1) = log10(5000); priors.Rmax(1) = log10(1000000);
    priors.Rmin(2) = log10(0.1); priors.Rmax(2) = log10(1000000);
else % num_layers = 1 
    priors.Rmin(1) = log10(0.1); priors.Rmax(1) = log10(1000000);
end

params.npt_init = 1; % initial number of floating nuclei
params.sigma_change_R = 4; 
params.sigma_move_depth = 10000;
params.sigma_birth_R = 6; % <<< Kept the increased value

% LESS IMPORTANT PARAMETERS (Define discretization grids)
params.thin = 1000; 
params.x_min = priors.layer_depths(1);
params.x_max = priors.depth_max;
x_temp = logspace(log10(params.x_min-1),log10(params.x_max),params.dis-1);
params.x = [0, x_temp]; % depth axis
params.y = linspace(min(priors.Rmin), max(priors.Rmax), params.dis); % R axis
params.frequency_x = linspace(log10(min(frequency)), log10(max(frequency)), params.dis);
% --- Define AR/Ph axes for plotting FORWARD MODELS ---
params.apparentRes_y = linspace(log10(min(dataAR_plot))-0.5, log10(max(dataAR_plot))+0.5, params.dis);
params.Phase_y = linspace((min(dataPh_plot)-20), (max(dataPh_plot)+20), params.dis);

%% %%%%%%%%%%%%%%%% RUN MULTIPLE CHAINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Starting %d MCMC chains...\n', num_chains);
results = cell(1, num_chains);

% Seed the *main* script once, outside the loop
rng('shuffle'); 
% Get a single "base" seed
baseSeed = randi(2^32 - 1);
fprintf('--- Using base seed: %d for chain generation ---\n', baseSeed);

parfor i = 1:num_chains
    % Give each worker a UNIQUE, deterministic seed
    chain_seed = baseSeed + i;
    fprintf('Running chain %d with seed %d...\n', i, chain_seed);
    % CALL run_single_chain_Zinvert function
    results{i} = run_single_chain_Zinvert(data, priors, params, i, chain_seed);
    fprintf('...Chain %d finished.\n', i);
end
fprintf('All chains finished. Combining results...\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINE AND PROCESS RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Combined_CI_density = results{1}.CI_density;
Combined_CI_density_FMAR = results{1}.CI_density_FMAR;
Combined_CI_density_FMPh = results{1}.CI_density_FMPh;
Combined_nnucleihist = results{1}.nnucleihist;
Combined_change_points = results{1}.change_points;
Combined_AV = results{1}.AV;
total_b = results{1}.b; 
best_like_overall = results{1}.like_best;
best_model_overall = results{1}.best_model_R;
Combined_cov = zeros(num_chains, params.nsample);
Combined_nnuclei = zeros(num_chains, params.nsample);
Combined_RMS = zeros(num_chains, params.nsample);
Combined_cov(1,:) = results{1}.cov;
Combined_nnuclei(1,:) = results{1}.nnuclei;
Combined_RMS(1,:) = results{1}.RMS;
for i = 2:num_chains
    Combined_CI_density = Combined_CI_density + results{i}.CI_density;
    Combined_CI_density_FMAR = Combined_CI_density_FMAR + results{i}.CI_density_FMAR;
    Combined_CI_density_FMPh = Combined_CI_density_FMPh + results{i}.CI_density_FMPh;
    Combined_nnucleihist = Combined_nnucleihist + results{i}.nnucleihist;
    Combined_AV = Combined_AV + results{i}.AV;
    total_b = total_b + results{i}.b;
    Combined_change_points = [Combined_change_points; results{i}.change_points];
    Combined_cov(i,:) = results{i}.cov;
    Combined_nnuclei(i,:) = results{i}.nnuclei;
    Combined_RMS(i,:) = results{i}.RMS;
    if results{i}.like_best < best_like_overall
        best_like_overall = results{i}.like_best;
        best_model_overall = results{i}.best_model_R;
    end
end
% Calculate final statistics
x = params.x; y = params.y;
R_edge = y(1:(params.dis-1)); depth_edge = x(1:(params.dis-1));
frequency_edge = params.frequency_x(1:(params.dis-1));
FM_edge = params.apparentRes_y(1:(params.dis-1)); % AR axis
FMPh_edge = params.Phase_y(1:(params.dis-1)); % Phase axis
AV = Combined_AV / total_b;
sizeCI = size(Combined_CI_density); len = sizeCI(1,1); wid = sizeCI(1,2);
CI_density_limitN = zeros(len,wid); total_sum = sum(Combined_CI_density, 2); 
for i = 1:len, if total_sum(i) > 0, CI_density_limitN(i,:) = Combined_CI_density(i,:) / total_sum(i); end, end
mode = zeros((params.dis-1),(params.dis-1));
for i=1:(params.dis-1), [~,I] = max(Combined_CI_density(i,:)); mode(i,I) = 1; end
R_mode = mode .* R_edge; R_mode_n0 = zeros((params.dis-1),1);
for i=1:(params.dis-1), idx = find(R_mode(i,:), 1); if ~isempty(idx), R_mode_n0(i,1) = R_mode(i, idx); else R_mode_n0(i,1) = NaN; end, end
R_mode_n0 = rmmissing(R_mode_n0); depth_mode = depth_edge(1:length(R_mode_n0)); 
sup = zeros(params.dis-1, 1); inf = zeros(params.dis-1, 1); cdf = zeros(wid, 1);
for i = 1:len, if total_sum(i) > 0, cdf = cumsum(Combined_CI_density(i,:)) / total_sum(i); inf_idx = find(cdf >= 0.025, 1, 'first'); sup_idx = find(cdf >= 0.975, 1, 'first'); if isempty(inf_idx), inf_idx = 1; end; if isempty(sup_idx), sup_idx = wid; end; inf(i) = R_edge(inf_idx); sup(i) = R_edge(sup_idx); else inf(i) = R_edge(1); sup(i) = R_edge(end); end, end
x2 = [x'; flipud(x')]; inbetween = [inf; flipud(sup)];

% --- Calculate 16% and 84% CIs (for 1-sigma equivalent) ---
sup_84 = zeros(params.dis-1, 1);
inf_16 = zeros(params.dis-1, 1);
cdf_16_84 = zeros(wid, 1);
for i = 1:len
    if total_sum(i) > 0
        cdf_16_84 = cumsum(Combined_CI_density(i,:)) / total_sum(i);
        inf_idx_16 = find(cdf_16_84 >= 0.16, 1, 'first');
        sup_idx_84 = find(cdf_16_84 >= 0.84, 1, 'first');
        
        if isempty(inf_idx_16), inf_idx_16 = 1; end
        if isempty(sup_idx_84), sup_idx_84 = wid; end
        
        inf_16(i) = R_edge(inf_idx_16);
        sup_84(i) = R_edge(sup_idx_84);
    else
        inf_16(i) = R_edge(1);
        sup_84(i) = R_edge(end);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT COMBINED RESULTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True R plot
R_edge_true = 10.^(R_edge); 
best_model_R_true = 10.^(best_model_overall);
inbetween_true = 10.^(inbetween); % This is the 95% CI
R_mode_n0_true = 10.^(R_mode_n0); 
AV_true = 10.^(AV);
% --- Convert 16/84 CIs to log10 space ---
inf_16_true = 10.^(inf_16);
sup_84_true = 10.^(sup_84);

figure('Name', 'Combined Resistivity PDF'); 
subplot(1,2,1);
contourf(R_edge_true, depth_edge, CI_density_limitN, 100, 'LineColor', 'none'); hold on; % Use CI_density_limitN
plot(R_mode_n0_true, depth_mode, 'k-', 'LineWidth', 2); 
plot(best_model_R_true, x, 'r-', 'LineWidth', 1.5); 
plot(inf_16_true, depth_edge, 'g--', 'LineWidth', 1.5); 
plot(sup_84_true, depth_edge, 'g--', 'LineWidth', 1.5); 
legend('Normalized PDF', 'Mode', 'Best Fit', '16/84% CI', 'Location', 'southwest'); 
title('Combined Resistivity PDF','FontSize',14);
ylabel('Depth (m)','FontSize',14); xlabel('R (Ohm.m)','FontSize',14); set(gca,'XScale','log');
set(gca,'XMinorTick','on'); set(gca,'Ydir','reverse'); colormap(cm_blue); colorbar; 
clim([0 max(CI_density_limitN(:))]); % Adjust color limit for normalized data

subplot(1,2,2); % Second plot (zoomed)
contourf(R_edge_true, depth_edge, CI_density_limitN, 100, 'LineColor', 'none'); hold on; % Use CI_density_limitN
plot(R_mode_n0_true, depth_mode, 'k-', 'LineWidth', 2); 
plot(best_model_R_true, x, 'r-', 'LineWidth', 1.5); 
plot(inf_16_true, depth_edge, 'g--', 'LineWidth', 1.5); 
plot(sup_84_true, depth_edge, 'g--', 'LineWidth', 1.5); 
legend('Normalized PDF', 'Mode', 'Best Fit', '16/84% CI', 'Location', 'southwest'); 
title('Combined Resistivity PDF','FontSize',14);
ylabel('Depth (m)','FontSize',14); xlabel('R (Ohm.m)','FontSize',14); set(gca,'XScale','log');
set(gca,'XMinorTick','on'); set(gca,'Ydir','reverse'); colormap(cm_blue); colorbar; 
clim([0 max(CI_density_limitN(:))]); % Adjust color limit
ylim([0 5000]);

%% FM PDF (log10 scale) - Plots AR/Ph

if data.running_mode == 1
    figure('Name', 'Combined Forward Model PDFs');
    subplot(1,2,1)
    imagesc(frequency_edge, FM_edge, Combined_CI_density_FMAR); hold on; % pdf plot
    % Plot data with errors (using dataAR_plot, rhoaerr_plot)
    y_plot_ar = log10(dataAR_plot);
    upper_bound_log = log10(dataAR_plot + rhoaerr_plot);
    lower_bound_log = log10(dataAR_plot - rhoaerr_plot);
    % Handle potential negative lower bounds before log10
    min_val_ar = min(dataAR_plot(dataAR_plot>0))/1e3; % small positive number
    lower_bound_linear_ar = dataAR_plot - rhoaerr_plot;
    lower_bound_linear_ar(lower_bound_linear_ar <= 0) = min_val_ar;
    lower_bound_log = log10(lower_bound_linear_ar);
    err_upper_delta_ar = upper_bound_log - y_plot_ar;
    err_lower_delta_ar = y_plot_ar - lower_bound_log;
    errorbar(data.frequency, y_plot_ar, err_lower_delta_ar, err_upper_delta_ar, 'k', 'markersize', 8, 'LineStyle', 'none');
    xlim([min(data.frequency), max(data.frequency)]);
    ylabel('Apparent Resistivity (Ohm.m)'); xlabel('Frequency (Hz)');
    set(gca,'XScale','log'); set(gca,'XDir','reverse'); set(gca,'XMinorTick','on'); 
    set(gca,'YDir','normal'); set(gca,'YMinorTick','on'); colormap(cm_blue); 
    colorbar; clim([0 10000]); ylim([1 5]); title('App. Res. PDF');

    subplot(1,2,2)
    imagesc(frequency_edge, FMPh_edge, Combined_CI_density_FMPh); hold on; % pdf plot
    errorbar( data.frequency, dataPh_plot, phierr_plot,'k','markersize',8, 'LineStyle', 'none'); 
    xlim([min(data.frequency), max(data.frequency)]);
    ylabel('Phase'); xlabel('Frequency (Hz)');
    set(gca,'XScale','log'); set(gca,'XDir','reverse'); set(gca,'XMinorTick','on'); 
    set(gca,'YDir','normal'); set(gca,'YMinorTick','on'); colormap(cm_blue); 
    colorbar; clim([0 10000]); title('Phase PDF');
end

%% %%%%%%%%%%%%% Plot Statistics %%%%%%%%%%%%%%%%
figure('Name', 'Chain Convergence'); subplot(2,1,1); hold on;
for i = 1:num_chains, semilogy(Combined_cov(i,:)); end
xlabel('iterations','FontSize',14); title('Data Misfit (All Chains)','FontSize',16); ylabel('Misfit');
subplot(2,1,2); hold on;
for i = 1:num_chains, plot(Combined_nnuclei(i,:)); end
line([params.burn_in params.burn_in],[0 priors.npt_max+priors.num_layers],'LineWidth',4,'Color',[1 0 0]);
xlabel('iterations','FontSize',14); title('Number of nuclei (All Chains)','FontSize',16);
ylabel('nNuclei'); ylim([0 priors.npt_max+priors.num_layers+5]);

%% %%%%%%%%%%%%% CONVERGENCE DIAGNOSTICS %%%%%%%%%%%%%%%%
% This section calculates and displays the Gelman-Rubin statistic (R_hat)
% and plots post-burn-in trace plots to check for mixing.
% 1. Calculate Gelman-Rubin Statistic (R_hat)
% We calculate it for two key parameters:
% - Misfit (log-likelihood)
% - Number of nuclei (a key model parameter)
if num_chains > 1
    % R_hat for Misfit (log-likelihood)
    % Using Combined_cov, which stores the likelihood ('like')
    R_hat_misfit = calculate_gelman_rubin(Combined_cov, params.burn_in);
    % R_hat for Number of Nuclei
    R_hat_nnuclei = calculate_gelman_rubin(Combined_nnuclei, params.burn_in);
    fprintf('Gelman-Rubin (R_hat) for N Nuclei: %.4f\n', R_hat_nnuclei); 
    if R_hat_misfit > 1.1 || R_hat_nnuclei > 1.1
        fprintf('WARNING: R_hat > 1.1. Chains may not have converged.\n');
    else
        fprintf('SUCCESS: R_hat < 1.1. Chains appear to have converged.\n');
    end
else
    fprintf('Cannot calculate Gelman-Rubin: only 1 chain was run.\n');
end
% 2. Plot Post-Burn-In Trace Plots (better for checking mixing)
post_burn_start = params.burn_in + 1;
post_burn_iter = post_burn_start:params.nsample;

figure('Name', 'Post-Burn-In Convergence');
subplot(2,1,1); hold on;
for i = 1:num_chains
    plot(post_burn_iter, Combined_cov(i, post_burn_start:end));
end
title('Post-Burn-In Trace: Misfit');
xlabel('Iteration'); ylabel('Misfit');
legend(arrayfun(@(i) sprintf('Chain %d', i), 1:num_chains, 'UniformOutput', false), 'Location', 'best');

subplot(2,1,2); hold on;
for i = 1:num_chains
    plot(post_burn_iter, Combined_nnuclei(i, post_burn_start:end));
end
title('Post-Burn-In Trace: Number of Nuclei');
xlabel('Iteration'); ylabel('nNuclei');
legend(arrayfun(@(i) sprintf('Chain %d', i), 1:num_chains, 'UniformOutput', false), 'Location', 'best');

%% %%%%%%%%%%%% Plot Marginals %%%%%%%%%%%%%%%%%%%
figure('Name', 'Combined Posteriors'); subplot(2,1,1);
hist(Combined_change_points, 500); title('Combined Probability of change points','FontSize',14);
xlabel('Depth (m)'); ylabel('Count');
subplot(2,1,2); bar(Combined_nnucleihist); 
title('Combined Posterior on number of nuclei ','FontSize',14);
xlabel('Number of Nuclei'); ylabel('Count');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SAVING WORKSPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('output.mat') %
fprintf('Combined results saved');