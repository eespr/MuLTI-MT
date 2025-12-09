%==========================================================================
% function [output] = run_single_chain_Zinvert(data, priors, params, chain_id)
% 
% MuLTI_MT
%
% This Matlab script runs a transdimensional MCMC inversion for 1-D MT data
% 
% Inverts using COMPLEX IMPEDANCE (Z_real, Z_imag)
% Calculates forward model histograms in AR/Ph space.
%
% The physical model consists of internal layers, each with an associated resistivity (R), sampled as log(R) defined by Voronoi nuclei.
% The domain is divided into a number of layers, num_layers (which could be one)
% each with its own prior distribution on R.
% Each layer has a special nuclei that cannot leave its layer ("confined" nuclei).
% npt is the number of "floating" nuclei that can change layer.
% It runs one MCMC chain and returns the results in a struct.
%
% Author: Siobhan Killingbeck, Swansea University (2025)
%==========================================================================

function [output] = run_single_chain_Zinvert(data, priors, params, chain_id, chain_seed)

%% %%%%%%%%%%%%%%%%%%%%%%% UNPACK INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack Z data and weights

% Data
frequency = data.frequency;
Z_data = data.Z_data; % Complex impedance data
weightingZ_real = data.weightingZ_real; % Real error/weight
weightingZ_imag = data.weightingZ_imag; % Imag error/weight
nd = data.nd;
running_mode = data.running_mode;

% Priors
num_layers = priors.num_layers;

% Parameters
burn_in = params.burn_in;
nsample = params.nsample;
dis = params.dis;
npt_init = params.npt_init;
sigma_change_R = params.sigma_change_R;
sigma_move_depth = params.sigma_move_depth;
sigma_birth_R = params.sigma_birth_R;
show = 1000; thin = params.thin;

% Discretization axes
x = params.x; y = params.y;
frequency_x = params.frequency_x;
apparentRes_y = params.apparentRes_y; % Axis for AR PDF plot
Phase_y = params.Phase_y;             % Axis for Phase PDF plot

% Seed RNG with the unique seed passed from the main script
rng(chain_seed); 
num=ceil((nsample-burn_in)*0.025/thin);  

%% %%%%%%%%%%%%%%%%%%%%%%% PRE-ALLOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=0; bb=0; AV=zeros(dis,1);
AB=0; AD=0; PB=0; PD=0; AcV=0; PV=0; AP=0; PP=0;
change_points=zeros(nsample*(priors.npt_max+num_layers),1);
cov=zeros(nsample,1); nnuclei=zeros(nsample,1);
CI_density=zeros((dis-1),(dis-1));
CI_density_FMAR = zeros((dis-1),length(frequency_x)-1); % AR PDF
CI_density_FMPh = zeros((dis-1),length(frequency_x)-1); % Phase PDF
nnucleihist=zeros(priors.npt_max+num_layers,1);
nuclei_depths=zeros(priors.npt_max+num_layers,1);
nuclei_R = zeros(priors.npt_max+num_layers,1);
RMS = zeros(nsample,1); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npt=npt_init;
for i=1:npt+num_layers, if i <= num_layers, if i == 1, top_of_layer = priors.depth_min; if num_layers > 1, bottom_of_layer = priors.layer_depths(1); else bottom_of_layer = priors.depth_max; end; elseif i < num_layers, top_of_layer = priors.layer_depths(i-1); bottom_of_layer = priors.layer_depths(i); else top_of_layer = priors.layer_depths(i-1); bottom_of_layer = priors.depth_max; end; nuclei_depths(i)= (bottom_of_layer + top_of_layer) / 2; else nuclei_depths(i)=priors.depth_min+rand*(priors.depth_max-priors.depth_min); end; layer = num_layers; for j = 1:num_layers - 1, if nuclei_depths(i) <= priors.layer_depths(j), layer = j; break; end; end; nuclei_R(i)=priors.Rmin(layer)+rand*(priors.Rmax(layer)-priors.Rmin(layer)); end 
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE INITIAL MISFIT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
like=0;
[thickness, R, priors_OK] = thicknesses_and_priors(nuclei_depths, nuclei_R, npt, num_layers, priors); 
R_log10 = R; % Keep log10 version for proposals if needed
R = 10.^(R); % Convert to linear for MT1D
if priors_OK == 1
    thickness = thickness(1,1:(length(thickness)-1)); 
    if running_mode == 1
      % MT1D returns AR, Phase, and Z
      [~, ~, Z_model] = MT1D(R,thickness,frequency); % Get complex Z 
      % --- Misfit Calculation using Z ---
      misfit_real = NaN(length(Z_model),1);
      misfit_imag = NaN(length(Z_model),1);
      for i = 1:length(frequency) 
          misfit_real(i,1) = real(Z_data(i,1)) - real(Z_model(i,1));
          misfit_imag(i,1) = imag(Z_data(i,1)) - imag(Z_model(i,1));
      end 
      like_real = nansum( (misfit_real).^2 ./(2 * weightingZ_real.^2) );
      like_imag = nansum( (misfit_imag).^2 ./(2 * weightingZ_imag.^2) );
      like = like_real + like_imag; % SUM of likelihoods

      % RMS Misfit (still based on Z)
      all_misfits = [misfit_real; misfit_imag]; 
      rms_total_norm = sqrt(mean(all_misfits .^ 2, 'omitnan')); 
    else
      like = 1; rms_total_norm = 1e99;
    end 
else % priors_OK == 0
    like = 1e99; % Assign high misfit if priors not met
    rms_total_norm = 1e99;
    Z_model = nan(size(frequency)); % Ensure Z_model exists but is NaN
end

like_best=1e99;
like_init=like;
rms = rms_total_norm;

%% %%%%%%%%%%%%%%%% START RJ-MCMC SAMPLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Chain %d: Total samples %i. Starting MCMC loop...\n', chain_id, nsample);

for s=1:nsample
    out=1;
    % Print statistics (Identical)
    if (mod(s,show)==0), number_of_samples = s; number_of_nuclei = npt; if (s>burn_in), fprintf('Chain %d | %7i     %5.2f         %5.2f       %5.2f     %5.2f\n', chain_id, s, 100*AcV/PV, 100*AP/PP, 100*AB/PB,100*AD/PD); end, end
     
    birth=0; move=0; death=0;
    nuclei_R_prop=nuclei_R; nuclei_depths_prop = nuclei_depths;
    like_prop = like;
    
    % --- Proposal Logic (Birth/Death/Move/Change R) ---
    % It proposes changes to nuclei_depths_prop and nuclei_R_prop
    % It calculates Prod_delta_prior, prob, move_prob
    if (mod(s,2)==0), if (s>burn_in), PV=PV+1; end; npt_prop = npt; ind=ceil(rand*(npt+num_layers)); nuclei_R_prop(ind) = nuclei_R(ind) + randn * sigma_change_R;
    else u=rand; 
        if (u<0.333), birth=1; if (s>burn_in), PB=PB+1; end; npt_prop = npt+1; nuclei_depths_prop(1:npt+num_layers) = nuclei_depths(1:npt+num_layers); nuclei_depths_prop(npt+num_layers+1) = priors.depth_min+rand*(priors.depth_max-priors.depth_min); 
            % --- Using Birth From Prior ---
            layer = num_layers; for j = 1:num_layers - 1, if nuclei_depths_prop(npt+num_layers+1) <= priors.layer_depths(j), layer = j; break; end, end; 
            nuclei_R_prop(npt+num_layers+1) = priors.Rmin(layer) + rand*(priors.Rmax(layer)-priors.Rmin(layer)); % Draw R from prior
            Prod_delta_prior = priors.Rmax(layer)-priors.Rmin(layer);
            prob = 1 / Prod_delta_prior; % Simplified Hastings ratio for uniform proposal
        elseif (u<0.666), death=1; if (s>burn_in), PD=PD+1; end; npt_prop = npt-1; ind=ceil(rand*npt)+num_layers; if npt > 0, nuclei_depths_prop(1:num_layers+npt-1) = [nuclei_depths(1:ind-1) ; nuclei_depths(ind+1:num_layers+npt)]; nuclei_R_prop(1:num_layers + npt-1)= [nuclei_R(1:ind-1) ; nuclei_R(ind+1:num_layers+npt)]; death_pt_R = nuclei_R(ind); death_pt_depth = nuclei_depths(ind); layer = num_layers; for j = 1:num_layers - 1, if death_pt_depth <= priors.layer_depths(j), layer = j; break; end, end; Prod_delta_prior = priors.Rmax(layer)-priors.Rmin(layer); prob = 1/Prod_delta_prior; else out=0; prob=1; Prod_delta_prior=1; end; % Cannot kill if only fixed nuclei left
        else move=1; if (s>burn_in), PP=PP+1; end; npt_prop = npt; ind=ceil(rand*(npt+num_layers)); if num_layers == 1 || ind > num_layers, nuclei_depths_prop(ind) = nuclei_depths(ind)+randn*sigma_move_depth; else if ind == 1, top_of_layer = priors.depth_min; bottom_of_layer = priors.layer_depths(1); elseif ind < num_layers, top_of_layer = priors.layer_depths(ind-1); bottom_of_layer = priors.layer_depths(ind); else top_of_layer = priors.layer_depths(ind-1); bottom_of_layer = priors.depth_max; end; nuclei_depths_prop(ind) = (bottom_of_layer-top_of_layer).*rand(1) + top_of_layer; end; 
            % Simplified move_prob, assumes symmetry or neglects it
            move_prob = 1; 
        end
    end
    %----------------------------------------------------------------------    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE MISFIT OF THE PROPOSED MODEL
    % Calculate misfit using Z_real and Z_imag
    if out==1
        like_prop=0;
        [thickness, R, priors_OK] = thicknesses_and_priors(nuclei_depths_prop, nuclei_R_prop, npt_prop, num_layers, priors);
        R_log10_prop = R; % Keep log10 version if needed later
        R = 10.^(R); % Convert for MT1D
        if priors_OK == 0 
            out = 0; like_prop = 1e99; % High misfit if priors fail
            Z_model_prop = nan(size(frequency)); % Ensure exists
        else
           thickness = thickness(1,1:(length(thickness)-1)); 
            if running_mode == 1
                [~, ~, Z_model_prop] = MT1D(R,thickness,frequency); % Get complex Z
                % --- Misfit Calculation using Z ---
                misfit_real = NaN(length(Z_model_prop),1);
                misfit_imag = NaN(length(Z_model_prop),1);
                for i = 1:length(frequency) 
                    misfit_real(i,1) = real(Z_data(i,1)) - real(Z_model_prop(i,1));
                    misfit_imag(i,1) = imag(Z_data(i,1)) - imag(Z_model_prop(i,1));
                end 
                like_real_prop = nansum( (misfit_real).^2 ./(2 * weightingZ_real.^2) );
                like_imag_prop = nansum( (misfit_imag).^2 ./(2 * weightingZ_imag.^2) );
                like_prop = like_real_prop + like_imag_prop; % SUM

                all_misfits = [misfit_real; misfit_imag]; 
                rms_total_norm = sqrt(mean(all_misfits .^ 2, 'omitnan'));
            else
                like_prop = 1; rms_total_norm = 1e99;
                Z_model_prop = nan(size(frequency));
            end
        end
    else % out was already 0 from proposal stage
        like_prop = 1e99; % Assign high misfit
        rms_total_norm = 1e99;
        Z_model_prop = nan(size(frequency));
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SEE WHETHER MODEL IS ACCEPTED %%%
    accept=0;
    if out == 0, Prod_delta_prior = 1; prob = 1; move_prob = 1; end % Avoid errors
    if (birth==1), if (rand<((1/(Prod_delta_prior*prob))*exp(log(out)-like_prop+like))), accept=1; if (s>burn_in), AB=AB+1; end; end
    elseif (death==1), if (rand<(Prod_delta_prior*prob*exp(log(out)-like_prop+like))), accept=1; if (s>burn_in), AD=AD+1; end; end
    elseif (move == 1), if (rand<(move_prob * exp(log(out)-like_prop+like))), accept=1; if (s>burn_in), AP=AP+1; end; end
    else if (rand<exp(log(out)-like_prop+like)), accept=1; if (s>burn_in), AcV=AcV+1; end; end; end
    
    % If accept, update the values (Store Z_model)
    current_Z_model = Z_model; % Store the model corresponding to 'like'
    if (accept==1)
        npt=npt_prop;
        nuclei_depths = nuclei_depths_prop;
        nuclei_R = nuclei_R_prop;
        like=like_prop;
        rms = rms_total_norm;
        current_Z_model = Z_model_prop; % Update Z_model to accepted one
    end
    
    % Calculate histograms AFTER acceptance
    % Use the current state (accepted or not) to calculate histograms
    if (s>burn_in)
        % Calculate R vs Depth PDF
        hist_R_temp=zeros(dis,1);
        for i=1:dis
            ind=whichnuclei(nuclei_depths(1:npt+num_layers),x(i),num_layers,priors);
            hist_R_temp(i)=nuclei_R(ind);
        end
        [N_temp]=histcounts2(x,hist_R_temp',x,y); 
        
        % --- Calculate AR/Ph PDFs from current_Z_model ---
        % Convert the current model's Z back to AR and Ph for histograms
        mua = 4*pi*10^-7;
        w = 2*pi*frequency;
        current_AR = (1./(w*mua)).*abs(current_Z_model).^2; 
        current_Ph = (atan2(imag(current_Z_model),real(current_Z_model))*(180/pi));

        % Interpolate onto histogram axes
        fmp_interpAR = interp1(log10(frequency),log10(current_AR),frequency_x,'linear','extrap');
        fmp_interpPh = interp1(log10(frequency),current_Ph,frequency_x,'linear','extrap');
        freq_no = linspace(1,dis,dis); % Dummy axis for histcounts2
        
        % Create the histograms using the AR/Ph axes
        [FM_AR_temp]=histcounts2(fmp_interpAR,freq_no,apparentRes_y,freq_no); 
        [FM_Ph_temp]=histcounts2(fmp_interpPh,freq_no,Phase_y,freq_no);       
        
        % Get edges
        R_edge=y(1:(dis-1)); depth_edge=x(1:(dis-1));
        FM_edge=apparentRes_y(1:(dis-1)); FMPh_edge=Phase_y(1:(dis-1));
        frequency_edge=frequency_x(1:(dis-1));
           
        % --- Sum the histograms ---
        CI_density = CI_density + N_temp; 
        CI_density_FMAR = CI_density_FMAR + FM_AR_temp; 
        CI_density_FMPh = CI_density_FMPh + FM_Ph_temp;
    end % if (s > burn_in) - histogram calculation
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collect samples
    if (s>burn_in)
        if (mod(s,thin)==0)
            b=b+1;
            for i=1:dis, ind=whichnuclei(nuclei_depths,x(i), num_layers,priors); AV(i,1)=AV(i,1)+nuclei_R(ind); end
            nnucleihist(npt+num_layers,1)=nnucleihist(npt+num_layers,1)+1;
            nuclei=nuclei_depths(1:npt+num_layers); nuclei=sort(nuclei);
            for i = 1:npt-1+num_layers, bb=bb+1; cp= (nuclei(i+1)+nuclei(i))/2; change_points(bb)=cp; end
        end
    end %if burn-in - sample collection
    
    % Store history
    cov(s)=like; 
    nnuclei(s)=npt+num_layers; 
    RMS(s) = rms; 
            
    % Get the best model (using 'like', which is based on Z misfit)
    if priors_OK == 1 && (like < like_best) % Check current model validity
        depths_best = nuclei_depths;
        R_best = nuclei_R;
        npt_best = npt;
        like_best = like;
    end   
    
end % the Sampling of the mcmc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating statistics of THIS chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_model_R = zeros(1,dis);
for i=1:dis
    ind=whichnuclei(depths_best(1:num_layers+npt_best),x(i),num_layers, priors);
    best_model_R(i)=R_best(ind);
end  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PACKAGE OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Chain %d: Packaging results.\n', chain_id);
output.CI_density = CI_density; output.CI_density_FMAR = CI_density_FMAR;
output.CI_density_FMPh = CI_density_FMPh; output.nnucleihist = nnucleihist;
output.change_points = change_points(1:bb); output.AV = AV; output.b = b;   
output.cov = cov; output.nnuclei = nnuclei; output.RMS = RMS;
output.best_model_R = best_model_R; output.like_best = like_best;

end % This 'end' closes the function