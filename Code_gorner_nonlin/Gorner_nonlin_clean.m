clear all; close all;

% Script last modified on 26.03.2020

% Output:
% edot_esti         posterior erosion rate    
% sigma_xy          map of posterior covariance
% red_sigma_xy      reduced sigma : sigma_post./sigma_prior
% res_xy            resolution for a given location
% res_dot           location of dot for which resolution is computed
% spread_xy         spread   
% e_p_e_pr_xy       difference posterior minus prediction with u_s : edot_esti - e_prior_us;


% Table of contents:
% O. Generalities: names and flags
% 0.1. Names of dates etc.
% O.2. Flags

% 1. Loading Gorner data
% 1.1. Load maps
% 1.2. Load source data
% 1.3. Load detrital data
% 1.4. Litho index to sample number

% 2. Parameters of the inverse model
% 2.1. Inverse parameters
% 2.2. Dimensions           /!\ no modifications from here on /!\
% 2.3. Total erosion 
% 2.4. Grid spacing etc.

% 3. Create the tracer matrix 

% 4. Inverse model 
% 4.1. Initialization 
% 4.2. Prior erosion rate
% 4.3. Model covariance matrix calculation 

% 5. Start iteration over different dates 
% 5.1. Detrital sample_i and gamma_obs 
% 5.2. Data covariance sample_i 
% 5.3. Start nonlinear iterations (stop when date misfit increases again) 
% 5.3.1.First step: start from prior or solve linear problem (flag) ---> solve, plot and update
% 5.3.2. Compute G ( G = d(g(m)) / dm ) 
% 5.3.3 Solve the non linear problem ---> solve, plot, update and compute misfit
% 5.4. Posterior uncertainties 
% 5.5. Save everything

% 6. Compute mean and median for dates without activation

%% This can all be modified

% 0. Generalities: names and flags
% 0.1. Names of data
load datenames
% datenames = {'24 July' '25 July' '27 July' '28 July' '30 July' ...
%     '1 Aug' '2_1 Aug' '2_2 Aug' '3 Aug' '4_1 Aug' '4_2 Aug'...
%     '5 Aug' '6 Aug' '7 Aug' '10 Aug' '11 Aug' '13 Aug' '14 Aug'...
%     '15 Aug' '27_1 Aug' '27_2 Aug' '27 3 Aug' '28 Aug' ' 29 Aug'...
%     '5 Sept' '7 Sept' '10 Sept' '11 Sept' '12 Sept' '13 Sept' '17 Sept'};
savename = 'Nonlinear_uniformprior';

%ind_GG = [19,59,71,18,60,61,62,63,68,69,70]; % gorner grenz samples
ind_GG = 1;

% 0.2 Flags
% 0.2.1. All dates or not
flag_alldates     = 1;  % if == 1: for all dates an erosion map will be computed

% 0.2.2. Uniform prior or u_s prior
flag_uniformprior = 1;  % if == 1: uniform prior, else: prior from U_s

% 0.2.3. First guess is prior or first guess is linear result (vals <0 = 0.01)
flag_firstprior   = 1;  % if == 1: start iterating from prior; else: start with solution linear model

% 0.2.4. Save everything
flag_saveall      = 0;  % if == 1: save all output

% 0.2.5. Sigma data from data measurements or self set parameter
flag_sigma_d_meas = 1;

% 0.2.6. Remove tracer data: specify tracer numbers to remove
flag_remove_data = 0;
to_remove = [6 7 18 21 23] - 2;

% 0.2.7. Plot mean result without activations: specify date numbers of act
flag_noactivations = 0;
datesnoactivation = [1 2 5 7 10 13 14 15 16 17 18 19 20 22 24 25 30]; % dates no activation

% 0.2.8. Plot data plots and misfit evolution
flag_plot_d_and_misfit = 0;

%1. Loading Gorner data
%1.1. Load maps
%[B,R] = geotiffread('litho_lowres.tif'); % Important: maps should be in same projection system, same extent and same cellsize
load('litholowres_watershed.mat');
[mask_GG,ggR] = geotiffread('mask_GGsamples.tif');
[DEM,R_DEM] = geotiffread('DEM_lowres.tif');
[vel_obs,R_v] = geotiffread('USurf_lowres.tif');
load geology_vector_shape

vel_obs(B<=0) = 0;
mask = vel_obs; % mask with glacier pixels
mask(B>0) = 1;

vel_obs(DEM>3600) = 0;  % no sliding above this altitude: cold ice
vel_obs = 0.5 .* vel_obs; % m/y
e_prior_us = 2.7 * 10^(-4) .* vel_obs.^2.02; % mm/y

%1.2. Load source data
load('source_B.mat')
% source = source_mean(:,:);  % cols: sample nums, rows: element nums (skip 1st 2 rows: silver)
source = source_percentage;
source(source==0) = 0.00001;    % get rid of 0 -> non linear model works with log() 
source(source==1) = 0.99;       
n_tracers = size(source,1);

%1.3. Load detrital data
load('DataSamples_B.mat');

%load('dataGunther_turbscaled.mat');
%datenames = dates;

Q_s = 100e6;                    % sediment load (kg) 2018
rho = 2650;                     % sediment density (kg/m3)
V_s = Q_s * (1/rho);            % total sediment volume (m3)

% 1.4. Litho index to sample index
S = unique(B);  % list of sample indices
S(1) = [];      % remove zero from list

B(B==6) = 4; % the lithology #6 is put as if it was #4 (expanding monte rosa granites)

% B((B==1)) = 60; %sample GR07
% B((B==2)) = 50; %sample GS08
% B((B==3)) = 40; %sample GS06
% B((B==4)) = 20; %sample GS02, GR02, GS03, GR03
% B((B==5)) = 30; %sample GR04, GS05
% B((B==6)) = 70; %sample GRmica
% B((B==7)) = 10; %sample GS01
% B = B/10;

% 2. Parameters of the model
% 2.1. Inverse parameters
sigma_cov = 3;  % standarddeviation model covariance (std:3)
L = 1200;       % critical distance covariance matrix (std 2500, best:1200)
d_sigma_val = 1e-2;

n_m = 20;                      % number of non-lin SD iterations
misfit = zeros(1,n_m);         % vector to store misfit
d_misfit = zeros(1,n_m);
m_misfit = zeros(1,n_m);
mu  = 1e-5;

%% Do not modify beyond here
% 2.2. Dimensions
[ny,nx] = size(B); % nx= number of cells in x-direction (#columns)
nn = nx*ny;        % number of gridcells of map
dx = R.CellExtentInWorldX;  % resolution in x direction (column)
dy = R.CellExtentInWorldY;
Ly = dy*(ny-1);
Lx = dx*(nx-1);
dA = dx * dy; % cell size (m2)

x2 = 0:dx:Lx;
y2 = Ly:-dy:0;
[X2,Y2] = meshgrid(x2,y2');

% 2.2. Initialise vectors
edot_v  = zeros(nn,1);  % erosion model to be filled in (vectorial form)
edot    = zeros(ny,nx); % erosion model to be filled in ('map' form)

obs_v           = reshape(mask,[nn,1]);
GG_v            = reshape(mask_GG,[nn,1]); GG_v(obs_v==0) = 0;
geol_v          = reshape(B,[nn,1]);

geol_v(obs_v==0)= [];

GG_v_inmask     = GG_v(obs_v == 1);
cols_GG         = find(GG_v_inmask == 1); % columns to select in A matrix

mask_size       = length(geol_v);
mask_size_GG    = length(cols_GG);

% 2.3. Total erosion
etot = (V_s / dA) * 1e3; % mm
%etot = 0.24 * mask_size; % mm
    
% 2.4. Grid spacing etc.
index=0;
x = zeros(1,nx*ny);
y = zeros(1,nx*ny);
for i = 1:nx
    for j = 1:ny
        index = index+1;
        x(1,index)= dx*(j-1);
        y(1,index) = dy*(i-1);
    end
end
x(obs_v==0) = [];
y(obs_v==0) = [];

% 3. Create the tracer matrix
A_0       = zeros(n_tracers, mask_size);    % initialization tracer matrix
for tr = 1:n_tracers
    A_v = geol_v;
    for geol_i = 1:length(S)
        A_v(A_v == S(geol_i)) = source(tr,geol_i);
    end
    A_0(tr,:) = A_v;
end
A_GG = A_0(:,cols_GG);

% 4. Inverse model
% 4.1. Initialization
if flag_alldates == 1
    ndates = length(datenames);
else
    prompt = 'Indicate the number of dates you want to evaluate';
    ndates = input(prompt);
end

solution_matrix = zeros(mask_size,ndates); % Empty matrix to store all solutions

mean_geol   = zeros(length(S),ndates); % Empty matrix to store mean erosion rates for different geologies
min_geol    = zeros(length(S),ndates);
max_geol    = zeros(length(S),ndates);
med_geol    = zeros(length(S),ndates);

% 4.2. Prior erosion rate
if flag_uniformprior == 1
    edot_prior_0 =  ones(mask_size,1) .* etot ./ mask_size;
else
    edot_prior_0 = reshape(e_prior_us,[nn,1]);
    edot_prior_0(obs_v==0) = [];
    edot_prior_0(edot_prior_0<0) = 0;
    edot_prior_0(edot_prior_0==0) = 0.00001;
end 

eps_prior_0 = log(edot_prior_0);
eps_prior_GG = eps_prior_0(cols_GG);


% 4.3. Model covariance matrix calculation
dist_0 = zeros(mask_size,mask_size);
cov_0  = zeros(mask_size,mask_size);
for i = 1:mask_size
    for j = i:mask_size
        dist_0(i,j)   = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        cov_0(i,j)    = sigma_cov^2*exp(-(dist_0(i,j)/L)^2);
        cov_0(j,i)    = cov_0(i,j); %model covariance
    end
end
dist_GG = dist_0(cols_GG,cols_GG);
cov_GG  = cov_0(cols_GG,cols_GG);

% 4.4 Set plotting pars
nS   = sqrt(ndates);
nCol = ceil(nS);        % # cols of subplots
nRow = nCol - (nCol * nCol - ndates > nCol - 1); % # rows of subplots

% 5. Start iteration over different dates
for data_i = 1:ndates
    disp(datenames{data_i})
    [r,c] = size(B);
    
    if any(ind_GG == data_i)
        eps_prior   = eps_prior_GG;
        A           = A_GG;
        cov         = cov_GG;
        dist        = dist_GG;
        n_pixels    = mask_size_GG;
        mask_v      = GG_v;
    else
        eps_prior   = eps_prior_0;
        A           = A_0;
        cov         = cov_0;
        dist        = dist_0;
        n_pixels    = mask_size;
        mask_v      = obs_v;
    end
    
    % 5.1. Detrital sample_i and gamma_obs
    d           = data_percentage(:,data_i);
    d           = d * etot; % to remove etot factor
    d(isnan(d)) = 0;
    d(d==0)     = 0.00001; % avoid problems with log()
    d(d==1)     = 0.99;
    
    gamma_obs   = log(d);
    
    % 5.2. Data covariance sample_i
    if flag_sigma_d_meas == 1
        d_sigma = data_sigma(:,1);
        d_sigma(d_sigma == 0) = 0.00001;
    else
        d_sigma = ones(size(data_percentage,1),1) * d_sigma_val;
%         d_sigma = d_sigma(3:end,1);
    end
    
    if flag_remove_data == 1
        d(to_remove) = [];
        gamma_obs(to_remove) = [];
        d_sigma(to_remove) = [];
    end
    
    if flag_remove_data == 1 && data_i == 1
        A(to_remove,:) = [];
        n_tracers = n_tracers - length(to_remove);
    end
    
    sigma_d             = d .* d_sigma;
    cov_d               = diag(sigma_d.^2);
    log_cov_d           = zeros(size(cov_d));
    log_cov_d(cov_d>0)  = log(cov_d(cov_d>0));
    cov_d               = log_cov_d; %data covariance
    
    
    % 5.3. Start nonlinear iterations (stop when date misfit increases again)
    for m = 1:n_m
        
        disp(['m is ' num2str(m)])
        
        % 5.3.1.First step: start from prior or solve linear problem (flag)
        if m == 1
            % Solve
            if flag_firstprior == 1              
                edot_post = exp(eps_prior); % Start from prior
            else
                H       = cov*A'*pinv((A*cov)*A'+cov_d);
                edot_post = exp(eps_prior)+H*(d-A*exp(eps_prior)); % Solve linear problem as initial guess
                edot_post(edot_post <= 0) = 0.0001;
            end
            
            % Update
            eps_m   = log(edot_post);
            gamma_m = log( A * edot_post );
            
            % Plot
             figure(1)
             e_plot = v_to_xy(edot_post,mask_v,ny,nx,1e-5);
             cmap = parula(200);
             cmap(1,:) = [1 1 1];
             imshow(e_plot,'Colormap',cmap,'DisplayRange',[-0.01 max(max(exp(eps_prior)))]);
             caxis([-0.01 max(max(exp(eps_prior)))])
             colorbar
             title(['Result for ' char(datenames(data_i)) ' after ' num2str(0) ' nonlinear iterations']);
             drawnow
             
        end
       
        
        % 5.3.2. Compute G ( G = d(g(m)) / dm )
        G = zeros(size(A));
        for i_i = 1:n_tracers
            for j_j = 1:n_pixels
                G(i_i,j_j) = (1/d(i_i)) * exp(eps_m(j_j)) * A(i_i,j_j) * (1/etot);
            end
        end
        
        
        % 5.3.3 Solve the non linear problem
        eps_m1 = eps_m + mu * ( cov * G' * cov_d^(-1) * ...
                 (gamma_m - gamma_obs) + (eps_m - eps_prior) );
        edot_post = exp(eps_m1);
        
        % Update
        eps_m   = eps_m1;
        gamma_m = log( A * edot_post );
        
        % Compute differences for computation misfit
        diff_gamma = gamma_m - gamma_obs;
        diff_eps   = eps_m   - eps_prior;
        
        % Compute misfit (total -, data -, model misfit)
        misfit(1,m) = diff_gamma' * pinv(cov_d) * diff_gamma + ...
                 diff_eps' * pinv(cov) * diff_eps;
        d_misfit(1,m) = diff_gamma' * pinv(log_cov_d) * diff_gamma;     % data misfit: difference observed vs modeled data
        m_misfit(1,m) = diff_eps' * pinv(cov) * diff_eps;               % model misfit: difference edot post vs edot prior
             
        % Plot posterior guess
        figure(1);
        e_plot = v_to_xy(edot_post,mask_v,ny,nx,1e-10);
        imshow(e_plot,'Colormap',cmap,'DisplayRange',[-0.01 max(max(exp(eps_prior)))]);
        caxis([-0.01 max(max(exp(eps_prior)))])
        colorbar
        title(['Result for ' char(datenames(data_i)) ' after ' num2str(m) ' nonlinear iterations']);
        drawnow

        % Stop iterating if data misfit increases
        if m > 1 && d_misfit(1,m-1) < d_misfit(1,m)
            disp('misfit kleiner')
            break
        end
        
    end
    
    % Plot solutions
    figure(2);
    subplot(nRow,nCol,data_i)
    imshow(e_plot,'Colormap',cmap,'DisplayRange',[-0.01 max(max(exp(eps_prior)))]);
    caxis([-0.01 max(max(exp(eps_prior)))])

    h = colorbar;
    set(get(h,'label'),'string','mm/yr');
    title([char(datenames(data_i)) ', ' num2str(m) ' n-l it.']);
    hold on
    yyaxis right
    plot(gla_out_ps_rel(:,1),gla_out_ps_rel(:,2),'k');
    hold off
    ax = gca;
    ax.PlotBoxAspectRatio = [size(DEM,2)/size(DEM,1) 1 1];
    

%     figure(); % same but 1 figure per sample:
%     imshow(e_plot,'Colormap', cmap, 'InitialMagnification', 1000);
%     caxis([-0.01 max(max(exp(eps_prior)))])
%     h = colorbar;
%     set(get(h,'label'),'string','mm/yr');
%     title([char(datenames(data_i)) ', ' num2str(m) ' n-l it.']);
%     hold on
%     yyaxis right
%     plot(gla_out_ps_rel(:,1),gla_out_ps_rel(:,2),'k');
%     hold off
%     ax = gca;
%     ax.PlotBoxAspectRatio = [size(DEM,2)/size(DEM,1) 1 1];
    
    if flag_plot_d_and_misfit == 1
    % Plot data obs vs data modelled
     figure(3);
     subplot(nRow,nCol,data_i);
     plot(exp(gamma_obs),exp(gamma_m),'o');
     refline(1,0);
     xlabel('gamma_{obs}'); ylabel('gamma_{mod}');
     title([char(datenames(data_i)) ', ' num2str(m) ' n-l it.']);
    
    % Plot data obs vs data modelled (histogram plot)
     figure(4);
     subplot(nRow,nCol,data_i);
     histogram(gamma_obs - gamma_m, length(d))
     title([char(datenames(data_i)) ', ' num2str(m) ' n-l it.']);
     title('gamma_{obs} - gamma_{m}')
    
     figure();
     plot(misfit)
     title('Misfit evolution')
    
     figure(3);
     plot(d_misfit)
     title('Data misfit evolution')
     
    end
    
    % 5.4. Posterior uncertainties
    H       = cov*A'*pinv((A*cov)*A'+cov_d);
    C_M     = cov - H*A*cov;
    R       = H*A;
    res_dot = edot_post*0; res_dot(20) = 1;
    spread  = edot_post*0;
    
    for i = 1:length(edot_post)
        for j = 1:length(edot_post)
            spread(i) = spread(i)+(dist(i,j)*R(i,j).^2);
        end
    end
    
    sigma_post   = sqrt(diag(C_M));
    sigma_prior  = sqrt(diag(cov));

    % 5.5. Save everything
    if flag_saveall == 1
  
%         sigma_post_v = mask_v; sigma_post_v(mask_v>0) = sigma_post;
%         sigma_prior_v = mask_v; sigma_prior_v(mask_v>0) = sigma_prior;
%         res20_v = mask_v; res20_v(mask_v>0) = R(20,:);
%         res_dot_v = mask_v; res_dot_v(mask_v>0) = res_dot;
%         spread_v = mask_v; spread_v(mask_v>0) = spread;
%         edot_post_v = mask_v; edot_post_v(mask_v>0) = edot_post;
        
        
        sigma_post_v = map_in_map_vector(sigma_post,mask_v);
        sigma_prior_v = map_in_map_vector(sigma_prior,mask_v);
        sigma_xy = v_to_xy(sigma_post,mask_v,ny,nx,0);
        red_sigma_xy = reshape((sigma_post_v./sigma_prior_v),[ny,nx]);
        res_xy = v_to_xy(R(20,:),mask_v,ny,nx,0);
        res_dot = v_to_xy(res_dot,mask_v,ny,nx,0);
        spread_xy = v_to_xy(spread,mask_v,ny,nx,0);
        edot_esti = v_to_xy(edot_post,mask_v,ny,nx,0); edot_esti(B==0)=NaN;
        e_p_e_pr_xy = edot_esti - e_prior_us;
        e_p_e_pr_xy(B==0)=0;
        
%         sigma_xy     = reshape((sigma_post_v),[ny,nx]);
%         red_sigma_xy = reshape((sigma_post_v./sigma_prior_v),[ny,nx]);
%         res_xy       = reshape(res20_v,[ny,nx]);
%         res_dot      = reshape(res_dot_v,[ny,nx]);
%         spread_xy    = reshape(spread_v,[ny,nx]);
%         edot_esti = reshape(edot_post_v,[ny,nx]);
%         edot_esti(B==0)=NaN;
%         e_p_e_pr_xy = edot_esti - e_prior_us;
%         e_p_e_pr_xy(B==0)=0;
%         e_max=max(max(edot_esti));
        
        savename_i = [savename '_' char(datenames(data_i)) '.mat'];
        save(savename_i, 'edot_esti', 'red_sigma_xy', 'sigma_xy','res_xy',...
            'res_dot','spread_xy','e_p_e_pr_xy');
    end
    
    % Save solution in solution matrix
    
    solution = map_in_obs_mask(edot_post,mask_v,obs_v);
    solution_matrix(:,data_i) = solution;

    [mean_geol(:,data_i),min_geol(:,data_i),...
        max_geol(:,data_i),med_geol(:,data_i)] = ...
        stats_for_geology(solution,geol_v,S);

end

% 6. Compute mean and median for dates without activation
if flag_noactivations == 1
    solutions = solution_matrix(:,datesnoactivation);
    solution_mean = mean(solutions,2);
    figure();
    e_plot = v_to_xy(solution_mean,obs_v,ny,nx,1e-10);
    imshow(e_plot,'Colormap',cmap,'DisplayRange',[-0.01 max(max(exp(eps_prior)))]);
    caxis([-0.01 max(max(exp(eps_prior)))])

    colorbar
    title('Mean erosion rate map for all dates without activation')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting real geological map in the end:
figure()
imagesc(B); axis image; title('Geological map')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Functions
function v = map_in_obs_mask(map_masked_v,mask_vector,obs_mask)
    v = mask_vector;  % vector of small or large mask vector
    v(mask_vector>0) = map_masked_v; 
    v(obs_mask==0) = []; 
end

function vv = map_in_map_vector(map_masked_v,mask_vector)
    vv = mask_vector;
    vv(mask_vector>0) = map_masked_v; 
end

function map_xy = v_to_xy(map_masked_v,mask_vector,ny,nx,zero_val)
    map_v = mask_vector; 
    map_v(mask_vector>0) = map_masked_v;
    map_v(mask_vector == 0) = -9999;
    map_v(map_v==0) = zero_val;
    map_xy = reshape(map_v, [ny nx]);
end

function [a,b,c,d] = stats_for_geology(e_post,geol_map,geol_indices)
    a = zeros(1,length(geol_indices));
    b = zeros(1,length(geol_indices));
    c = zeros(1,length(geol_indices));
    d = zeros(1,length(geol_indices));
    
    for geol_ind = 1:length(geol_indices)
        geol = geol_indices(geol_ind);
        e_geol_i = e_post(geol_map == geol);
            if isempty(e_geol_i) == 0
                a(1,geol_ind) = mean(e_geol_i);
                b(1,geol_ind) = min(e_geol_i);
                c(1,geol_ind) = max(e_geol_i);
                d(1,geol_ind) = median(e_geol_i);
            end
    end
end
