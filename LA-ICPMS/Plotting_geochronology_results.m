%{

Bruno Belotti (bruno.belotti@unil.ch) - Master's Thesis 2020/2021 processing - UNIL
Last modified: 12.03.2021
Plotting geochronology results for single-grain U-Pb LA-ICPMS results.

This is just a self-written tool and should not be used as a reliable interpretation key for 
geochronology results!

%}
%% INSTRUCTIONS:
%{
If you want to use this code on my original data, just run it without doing anything.

If you want to use this code to have a look at your personal data, do as follows:

Format of input data:
For each sample (LAMTRACE output), create one '.mat' file:
Each file must contain one nX10 matrix, as follows:
                                                                                      (optional)
   1         2         3         4         5         6         7         8       9        10
206/238 | 206/238 | 207/235 | 207/235 | 207/235 | 207/235 | 206/238 | 206/238 | Rho | discordance
  Age   | ?2sigma |   Age   | ?2sigma |  Ratio  | 1sigma% |  Ratio  | 1sigma% |     |      % 
   .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |  .  |      .
   .         .         .         .         .         .         .         .       .         .
   .         .         .         .         .         .         .         .       .         .

To do so, copy-paste the correct columns from the LAMTRACE outputs directly into an 
empty matlab variable (open it in the workspace): 
MySampleName = [];

Then rename it and save it as a .mat file:
save('MySampleName','MySampleName')

Note that data is structured according to LAMTRACE's reduction outputs (Jackson, S.E., 2008).
'Discordance %' can be used to roughly estimate result quality and is calculated as follows:

|1 - ([206Pb/238U AGE] / [207Pb/235U AGE]) *100|

Full column names (for info only):
DataFormat = {'206Pb/238U Age' '206Pb/238U ?2sigma' '207Pb/235U Age' '207Pb/235U ?2sigma'... 
              '207Pb/235U Ratio' '207Pb/235U 1 sigma %' '206Pb/238U Ratio' '206Pb/238U 1 sigma %'...
              'Rho' 'discordance %'}; save('DataFormat','DataFormat');
%}

%% Defining parameters: --------------------------------------------------------------------
clear all; clc; close all;
% Once the data are ready, follow the next steps and fill in all of the following details (modifying
% the arrays' sizes when necessary):

% 1) In the current directory (of this '.m' file), create a folder called 'data' and put all your 
%    samples' '.mat' data matrixes (in the above specified format) in it.

% 2) OPTIONAL - Specify discrodia percentage threshold (!!! All higher results will be removed, 
%    put '100' if you want to keep all results !!!):
disc_perc = 100; % [%]

% 3) If you don't want to use a 'discordance %' column, leave the 10th column empty 
%    and replace this '1' with a '0':
discord_col = 0;

% 4) Specify maximum ages (in Ma) to be considered for each sample:
%    (if you want all ages [and if your samples are from planet earth], set '4600')
%       Case a) If same max age for all samples:
        maxages = ones(1,17) * 4600; % <-- modify array length according to your number of samples,
%                                          and modify maximum age.

%       Case b) If different values for each sample:
%     maxages = [4600 4600 4600 4600 ...
%                4600 4600 4600 4600 ...
%                4600 4600 4600 4600 ...
%                4600 4600 4600 4600 ...
%                4600 ]; % [Ma] <-- modify array length according to your number of samples!

% 5) Define age bin intervals to be displayed in plots for each sample:
%       Case a) If same age bins for all samples:
        agebins = ones(1,17) * 5; % <-- modify array length according to your number of samples,
%                                       and modify bin size.
%       Case b) If different age bins for each sample:
%         agebins = [5 5 5 5 ...
%                    5 5 5 5 ...
%                    5 5 5 5 ...
%                    5 5 5 5 ...
%                    5 ]; % [Ma] <-- modify array length according to your number of samples!

% 6) Decide whether to plot single points or error bars (not recommended for big datasets). 
%    '0' = no errorbars; '1' = errorbars.
errellips   = 1;
scalefactor = 1.7; % 2.44 for 95% confidence from sigma error ellipses.

% ----------------------------------------------------------------------------------------
%% Data selection:
files = dir('data/*.mat');
for i=1:length(files)
    eval(['load data/' files(i).name]);
    samplenames{i} = files(i).name(1:end-4);
end
sample = 1:length(samplenames);
samplenames = [samplenames;num2cell(sample)]; disp(samplenames);
% One number was assigned to each sample (for indexing purposes).

% Now, we create a cell array containing all samples data according to their assigned indexes:
for i = sample
    MySamples{i} = eval(samplenames{1,i});
end; clear i;

% Here, we check if the discordance percentage column must be ignored, if not, we use the
% discordance percentage threshold and maximum ages to remove highly discordant results
% and results that are too old:
for i = sample
    if discord_col == 0
        MySamples{i}(:,10) = 0;
    end
    ind_disc{i}  = MySamples{i}(:,10) < disc_perc; % index of each 'concordant' result
    data{i}      = MySamples{i}(ind_disc{i},:);    % data are the new matrixes w/o disc. results
    ind_2old{i}  = data{i}(:,1) < maxages(i);      % index of results with age < maxages
    data{i}      = data{i}(ind_2old{i},:);         % updating data matrixes
end; clearvars i ind_disc ind_2old;

% Now our data is "clean" and we can proceed to plotting.

%% Plotting:

for i = sample
    figure(i)
% Concordia diagrams:
    subplot(2,1,2)
  if errellips == 0 % version without errorbars:
    t{i} = 0:1:round(max(data{i}(:,1))); % this sets the axis limits according to the sample's 
                                         % age distribution.
    Pb206_U238{i} = exp((1.55125e-10*1e6)*t{i})-1;   % these create ratios along the correct 
    Pb207_U235{i} = exp((9.8485e-10*1e6)*t{i}) -1;   % time interval.
    plot(Pb207_U235{i},Pb206_U238{i},'k'); hold on;
    scatter(data{i}(:,5),data{i}(:,7),'x','r');
        % Finding age indexes:
        interv = round(max(data{i}(:,1))/6,-2) : ...
                 round(max(data{i}(:,1))/6,-2) : ...
                 round(max(data{i}(:,1)));
        for k = 1:length(interv)
            ind{k} = t{i} == interv(k);
            scatter(Pb207_U235{i}(ind{k}),Pb206_U238{i}(ind{k}),150,'.','k');
            text(Pb207_U235{i}(ind{k}),Pb206_U238{i}(ind{k})-(max(Pb206_U238{i})/20), ...
                 num2str(interv(k))+" Ma", 'Color',[0 0 0],'FontSize',7);
        end
    
    hold off; grid on;
    title(samplenames{1,i} + "- Concordia diagram");
    xlabel('^{207}Pb/^{235}U'); ylabel('^{206}Pb/^{238}U');
    
  else % version with error ellipses: ------------------------------------------------------------
    t{i} = 0:1:round(max(data{i}(:,1))); % this sets the axis limits according to the sample's 
                                         % age distribution.
    Pb206_U238{i} = exp((1.55125e-10*1e6)*t{i})-1;   % these create ratios along the correct 
    Pb207_U235{i} = exp((9.8485e-10*1e6)*t{i}) -1;   % time interval.
    
    subplot(2,1,2); hold on
    plot(Pb207_U235{i},Pb206_U238{i},'k'); hold on;
    
    aspectratio{i} = (max(data{i}(:,5)-min(data{i}(:,5)))/max(data{i}(:,7)-min(data{i}(:,7))));
    for u = 1:length(data{i}(:,1)) % repeat once for each ratio combination (in each sample) 
%         covmat{u} = cov(data{i}(:,5),data{i}(:,7)); % ratios covariance matrix.
        cova{i}(u)   = data{i}(u,9)*data{i}(u,6)*data{i}(u,8);
        covmat{i}{u} = [data{i}(u,6)^2 cova{i}(u); cova{i}(u) data{i}(u,8)^2];
        eigval{i}{u} = eig(covmat{i}{u}); % calculating lengths of ellipse axes (= square roots 
                                          % of the eigenvalues of the covariance matrix).
        axlen{i}{u}  = sqrt(eigval{i}{u})*scalefactor;

        % Finding ellipse angles:
        theta{i}{u} = (1/2)*atan((1/aspectratio{i})*((2*covmat{i}{u}(1,2))/(covmat{i}{u}(1,1) ...
                                                                           -covmat{i}{u}(2,2))));
        % Creating the ellipses:
        % Rotation matrix:
        R{i}{u}    = ([cos(theta{i}{u}) -sin(theta{i}{u}); sin(theta{i}{u}) cos(theta{i}{u})]);        

        xc{i}{u} = data{i}(u,5);
        yc{i}{u} = data{i}(u,7);
        axiswh{i}{u} = [data{i}(u,6)/100*data{i}(u,5)*2; data{i}(u,8)/100*data{i}(u,7)*2];

        a{i}{u}  = data{i}(u,6)/100*data{i}(u,5);
        b{i}{u}  = data{i}(u,8)/100*data{i}(u,7);
        m{i} = length(t{i});
        theta0{i} = linspace(0,2*pi,m{i});

        for k = 1:m{i}
            x{i}{u}(k) = a{i}{u} * cos(theta0{i}(k));
            y{i}{u}(k) = b{i}{u} * sin(theta0{i}(k));
        end
        rCoords{i}{u} = R{i}{u} * [x{i}{u}; y{i}{u}];
        xr{i}{u} = rCoords{i}{u}(1,:);
        yr{i}{u} = rCoords{i}{u}(2,:);

        plot(xr{i}{u}+xc{i}{u},yr{i}{u}+yc{i}{u},'r')
    end % clearvars covmat coveig theta axlen;
    
    
    %---------------------------------------------------------------------------------------------
        % Finding age indexes:
    interv = round(max(data{i}(:,1))/6,-2) : ...
             round(max(data{i}(:,1))/6,-2) : ...
             round(max(data{i}(:,1)));
    for k = 1:length(interv)
        ind{k} = t{i} == interv(k);
        scatter(Pb207_U235{i}(ind{k}),Pb206_U238{i}(ind{k}),150,'.','k');
        text(Pb207_U235{i}(ind{k}),Pb206_U238{i}(ind{k})-(max(Pb206_U238{i})/20), ...
             num2str(interv(k))+" Ma", 'Color',[0 0 0],'FontSize',7);
    end
    
        hold off; grid on;
        title(samplenames{1,i} + "- Concordia diagram");
        xlabel('^{207}Pb/^{235}U'); ylabel('^{206}Pb/^{238}U');  
    end
    
% Age signatures:
    subplot(2,1,1)
    [counts{i}, binCenters{i}] = hist(data{i}(:,1), round(max(data{i}(:,1))/agebins(i)));
    bar(binCenters{i}, counts{i}, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [0.5 0.5 0.5]); hold on;
    ylabel('Count'); xlabel('Age [Myr]'); 
    yyaxis right
    set(gca,'ytick',[]);
    % Choosing 206Pb/238U for ages <800 Ma, and 207Pb/235U for ages >800 Ma:
    agedata{i}(data{i}(:,1) <= 800) = data{i}(data{i}(:,1) <= 800);
    agedata{i}(data{i}(:,1)> 800) = data{i}(data{i}(:,1)> 800);
    [f{i},xi{i}] = ksdensity(agedata{i}); % using new data for plot (mix between 206/238 and 
                                          % 207/235 ages) according to resulting ages.
    plot(xi{i},f{i},'k');
    title(samplenames{1,i} + "- U-Pb age distribution");
    set(gca,{'ycolor'},{'k'});
    hold off; grid on;
end


%% BEDROCK Multi-plots along same time intervals:

bcolors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560] ...
           [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840]};

figure(3)
subplot(3,1,3)
for i = [1 2 3 4 5 6]
    agedata{i}(data{i}(:,1) <= 800) = data{i}(data{i}(:,1) <= 800);
    agedata{i}(data{i}(:,1)> 800) = data{i}(data{i}(:,1)> 800);
    
    [counts{i}, binCenters{i}] = hist(data{i}(:,1), round(max(data{i}(:,1))/agebins(i)));
    hold on;
    bar(binCenters{i}, counts{i}, 'FaceAlpha', 0.5,'FaceColor',bcolors{i}, ...
        'EdgeColor', [0.5 0.5 0.5]); legend(samplenames{1,[1 2 3 4 5 6]});
    axis([0 800 0 15]); grid on;
    xticks([0:50:800])
    title("Bedrock U-Pb age distribution");
    ylabel('Count'); xlabel('Age [Myr]'); 
    colormap jet 
    
end
hold off;

% WATER Multi-plots along same time intervals:
subplot(3,1,1)
for i = 11
    agedata{i}(data{i}(:,1) <= 800) = data{i}(data{i}(:,1) <= 800);
    agedata{i}(data{i}(:,1)> 800) = data{i}(data{i}(:,1)> 800);
    
    [counts{i}, binCenters{i}] = hist(data{i}(:,1), round(max(data{i}(:,1))/agebins(i)));
    hold on;
    bar(binCenters{i}, counts{i}, 'FaceAlpha', 0.5,'FaceColor',bcolors{1}, 'EdgeColor', [0.5 0.5 0.5]); 
    legend(samplenames{1,11});
    axis([0 800 0 30]); 
    grid on;
    title("Suspended U-Pb age distribution");
    ylabel('Count'); % xlabel('Age [Myr]'); 
    colormap jet
    set(gca,'xticklabel',{[]})
    
end
hold off;

% SAND Multi-plots along same time intervals:

subplot(3,1,2)
for i = 7
    agedata{i}(data{i}(:,1) <= 800) = data{i}(data{i}(:,1) <= 800);
    agedata{i}(data{i}(:,1)> 800) = data{i}(data{i}(:,1)> 800);
    
    [counts{i}, binCenters{i}] = hist(data{i}(:,1), round(max(data{i}(:,1))/agebins(i)));
    hold on;
    bar(binCenters{i}, counts{i}, 'FaceAlpha', 0.5,'FaceColor',bcolors{3}, 'EdgeColor', [0.5 0.5 0.5]); 
    legend(samplenames{1,7});
    axis([0 800 0 40]); 
    grid on;
    title("Deposited U-Pb age distribution");
    ylabel('Count'); % xlabel('Age [Myr]'); 
    colormap jet
    set(gca,'xticklabel',{[]})
    
end
hold off;

%% WATER plots for specific time intervals (GW27+GW29 = 05.06.2019 & other GWs = 28.06.2019):
figure()
subplot 211
for i = [16 17]
    agedata{i}(data{i}(:,1) <= 800) = data{i}(data{i}(:,1) <= 800);
    agedata{i}(data{i}(:,1) >  800) = data{i}(data{i}(:,1) >  800);
    
    [counts{i}, binCenters{i}] = hist(data{i}(:,1), round(max(data{i}(:,1))/agebins(i)));
    hold on;
    bar(binCenters{i}, counts{i}, 'FaceAlpha', 0.5,'FaceColor',bcolors{1}, 'EdgeColor', [0.5 0.5 0.5]); 
    legend(samplenames{1,16:17});
    axis([0 800 0 15]);
    grid on;
    title("Suspended U-Pb age distribution");
    ylabel('Count'); % xlabel('Age [Myr]'); 
    colormap jet
    
end

subplot 212
for i = 12:15
    agedata{i}(data{i}(:,1) <= 800) = data{i}(data{i}(:,1) <= 800);
    agedata{i}(data{i}(:,1)> 800) = data{i}(data{i}(:,1)> 800);
    
    [counts{i}, binCenters{i}] = hist(data{i}(:,1), round(max(data{i}(:,1))/agebins(i)));
    hold on;
    bar(binCenters{i}, counts{i}, 'FaceAlpha', 0.5,'FaceColor',bcolors{2}, 'EdgeColor', [0.5 0.5 0.5]); 
    legend(samplenames{1,12:15});
    axis([0 800 0 15]); 
    grid on;
    title("Suspended U-Pb age distribution");
    ylabel('Count'); % xlabel('Age [Myr]'); 
    colormap jet
    
end
hold off;

