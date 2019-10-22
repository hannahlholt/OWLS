% Photocurrent_EUVOP_v2.m is a program that can be used to investigate the
% occulatation performance of the OWLS EUVOP instrument. This includes making a model
% atm using MSIS and finding the extinction ratio for the specified
% channels. See the 'Plot Functions section below'. 

% VERSION 2.0           Compares MAVEN EUVM Al/Nb/C and A/Ti/C filters with
%                       with Zr on EUVOP
% VERSION 3.0           Adds ability to use AlMg filter to compute O and N2
%                       denisties separately

close all;


% EDIT THESE.....
%-------------------------------------------------------------------------
IDLfile_dir = '/Users/hannahholt/Documents/IDLfiles/';  %directory with .sav files
irradiance_file = 'FISM_EVE_daily_merged_earth_v8r1_01_2007_on.sav';    % FISM irradiances
% MAVEN filter response func. EUV_OP = Zr and AlMg, euvma = Al/Nb/C and euvmb = A/Ti/C
resp_func_file = {'EUV_OP_resp_function.sav', 'euvma_response_cal.sav', 'euvmb_response_cal.sav'};                    
cross_sec_file = 'photon_cross_sections.sav';

% used only as a check 
LYRA_Zr_file = 'N_col_o_2014-12-20_41.9064_km.sav';                     

% Plot Functions
% --------- SIMULATIONS ---------
plot_densities = 0;             % To plot MSIS mass density outputs
plot_irradiance = 0;            % To plot an example irradiance
plot_cross_section = 1;         % To plot the photon cross sections
plot_LYRA_Zr = 0;               % To plot the example data set from the LYRA Zr channel
plot_respFunc = 1;              % To plot the response functions for each channels
plot_photocurrent = 0;          % To plot the photocurrents
plot_col_densities = 0;         % To plot the integrate LOS column number densities
plot_extinction = 0;            % To plot the extinction ratio for the equinoxes

% --------- INVERSIONS ----------
plot_N_convergence_LYRA = 0;    % To plot the converged O+N2 LOS column number densities with LYRA actual data (Zr filter)
plot_inverted_dens_LYRA = 0;    % To plot the O+N2 LYRA computed number density and compare with MSIS output (Zr filter)
plot_N_convergence_Zr = 1;      % To plot the converged O+N2 LOS column number densities from MSIS
plot_inverted_dens_Zr = 1;      % To plot the O+N2 computed number density and compare with direct MSIS output
plot_N_convergence_AlMg = 1;    % To plot the converged O,N2 LOS column number densities from MSIS
plot_inverted_dens_AlMg = 1;    % To plot the O,N2 computed number density and compare with direct MSIS output

plotyear = 2010;                % which year do you want to plot irradiances
plotlambda = 10;                % which wavelength do you want to plot [nm]
which_date = 1;                 % which date you want to plot (1 = spring eqnx, 2 = fall eqnx)

R_sun = 6.957E8;                % Solar radius [m]
b = 1.49597870700E11;           % avg distance from Earth to Sun (m)
zeta_max = atan(R_sun/b);       % Ang from center of Sun to edge of disk (AS SEEN FROM EARTH!) [rad]
N_zeta_pts = 100;
dzeta = zeta_max/N_zeta_pts;
num_filters = length(resp_func_file) + 1;
%-------------------------------------------------------------------------

%% Use MSIS to create model spherical atmosphere
%--------------------------------------------------------------------------
yr = 2010;
spr_eqnx = [yr, 3, 10];          % 2010  spring equinox [yr, mn, dy]
fall_eqnx = [yr, 9, 22];
amu = 1.66054e-27;               % [kg/amu]
R = 6371.0E3;                    % volumetric radius of Earth [m]
massO = 16*amu;
massO2 = 32*amu;
massN2 = 28*amu;

minz_atm = 50E3;                % min alt at tangent point, below which OWLS cannot take occulation measurement [m] 
maxz_atm = 1000E3;              % max alt cutoff for atmosphere [m] 
z_sc = 500E3;                   % (circular) orbital altitude of OWLS [m] 
    
minz_occ = minz_atm;            % min alt at tangent point, below which OWLS cannot take occulation measurement [m]                             
maxz_occ = z_sc;                % max alt at tangent point, equal to z_sc because highest tangent ray point is when s/c is directly overhead   

dz = 1000;                      % altitude spacing [m]

z_atm = (minz_atm:dz:maxz_atm)';                % atmospheric altitude profile [m], dz spacing
z_occ = (minz_occ:dz:maxz_occ)';                % tangent point altitude values [m], dz spacing
num_z_atm = length(z_atm);
num_z_occ = length(z_occ);

lat = -82 * ones(num_z_atm,1);          % latitude [deg]                   (z x 1)
lon = -90 * ones(num_z_atm,1);          % longitude [deg]                  (z x 1)  
UT = 12 * ones(num_z_atm,1);            % we want to measure at dawn terminator = 90W of UTC noon  [s]
year = yr * ones(num_z_atm,1);

% day of year for spring and fall eqnx
doy1 =  jd2doy(cal2jd(spr_eqnx(1), spr_eqnx(2), spr_eqnx(3)));   
doy2 = jd2doy(cal2jd(fall_eqnx(1), fall_eqnx(2), fall_eqnx(3)));
doy = [doy1 .* ones(num_z_atm,1), doy2 .* ones(num_z_atm,1)];
num_dates = length(doy(1,:));

% define density outputs of MSIS-00
rho_tot = zeros(size(doy));                               %(rows = alt, cols = dates)
rho_O2 = rho_tot; rho_N2 = rho_tot; rho_O = rho_tot;

% Call MSIS and get densities at every altitude point and every day
% output of atmosnrlmsise00:
% T    :   (n x 2) The first column is exospheric temp, the second is temperature at altitude  [K]
% rho  :   (n x 9) densities of species, columns are as follows:
%          coln 1: He               [1/m^3]
%          coln 2: O                [1/m^3]
%          coln 3: N2               [1/m^3]
%          coln 4: O2               [1/m^3]
%          coln 5: Ar               [1/m^3]
%          coln 6: Total            [kg/m^3]        (sum of everything except O_anomalous)
%          coln 7: H                [1/m^3]
%          coln 8: N                [1/m^3]
%          coln 9: O_anomalous      [1/m^3]
                
for i=1:num_dates
    [T, rho] = atmosnrlmsise00(z_atm, lat, lon, year, doy(:,i), UT);      % call MSIS
    rho_O(:,i) = rho(:,2);         % number density of species [1/m^3] for each altitude [row], and date [column]
    rho_N2(:,i) = rho(:,3);
    rho_O2(:,i) =  rho(:,4);
    rho_tot(:,i) = rho(:,6);        % mass density [kg/m^3]     
end


% check the densities
if plot_densities
    mass_dens = figure;
    semilogx(massO.*rho_O(:,1), z_atm/1000, massN2.*rho_N2(:,1), z_atm/1000, massO2.*rho_O2(:,1), z_atm/1000, rho_tot(:,1), z_atm/1000, 'k', 'linewidth', 2);   % mass dens
%     semilogx(rho_O(:,1), z_atm/1000, rho_N2(:,1), z_atm/1000, rho_O2(:,1), z_atm/1000, 'linewidth', 2);       % num dens 
    legend('O', 'N2', 'O2', 'Total');
    grid on;
    xlabel('Mass Densities [kg/m^3]');
    ylim([50 500])
    ylabel('Altitude [km]');
    title('NRLMSIS-00 Mass Densities for 3/20/2010');
end

%% READ IN the Irradiance Data
%--------------------------------------------------------------------------
data = restore_idl([IDLfile_dir, irradiance_file]);

lambda = data.FISM_WV;              % wavelength of spectrum [nm]                    (nx1)
dates = data.DAY_AR;                % date , format [yyyydoy]                        (mx1)
irradiance = (data.FISM_PRED);      % FISM model irradiance predictions [W/m^2/nm]   (mxn)
                                    % for each day and each wavelength
                                    % rows = day, columns = irradiance for
                                    % specific lambda
num_lambda = length(lambda);        % number of different wavelengths for irradiance data                                        
                                    
% extract years and day of years
[year, doy] = extract_yearDOY(dates);
plotdoy = doy((year == plotyear)) ;
[dis, index] = min(abs(lambda - plotlambda));                               

if plot_irradiance
    irr = figure;
    subplot(2,1,1)
    plot(lambda, irradiance(1,:));
    xlabel('Wavelength [nm]');
    ylabel('FISM Predicted Irradiance [W/m^2/nm]');
    title(['FISM Predicted Irradiance at 1AU on ', num2str(dates(1))]);
    grid on;
    
    subplot(2,1,2)
    plot(plotdoy, irradiance((year == plotyear), index));
    xlabel('DOY');
    xlim([1 365])
    ylabel([num2str(lambda(index)), ' nm Irradiance [W/m^2/nm]']);
    title(['FISM Predicted Irradiance at 1AU for ', num2str(plotyear)]);
    grid on;
    
end

%% READ IN the photon cross sections
%--------------------------------------------------------------------------

data = restore_idl([IDLfile_dir, cross_sec_file]);

lambda_CS = data.PHOTO.ANGSTROMS ./10;      % wavelength for the cross sections [nm] 
sigma_O = data.PHOTO.O3P.XSECTION(1,:)/1E4;     % atomic oxygen cross section [cm^2]
sigma_N2 = data.PHOTO.N2.XSECTION(1,:)/1E4;     % molecular nitrogen cross section [cm^2]
sigma_O2 = data.PHOTO.O2.XSECTION(1,:)/1E4;     % moelcular oxygen cross section [cm^2]

% interpolate the cross sections to find at every lambda
% value from the irradiance data (since it is much coarser than the cross section wavelengths)
sigma_O_itp = interp1(lambda_CS, sigma_O, lambda);
sigma_N2_itp = interp1(lambda_CS, sigma_N2, lambda);
sigma_O2_itp = interp1(lambda_CS, sigma_O2, lambda);


if plot_cross_section
    cross_sec = figure;
    semilogy(lambda_CS, sigma_O, 'b', lambda_CS, sigma_N2, 'r', lambda_CS, sigma_O2, 'g');
    xlim([lambda_CS(1), lambda_CS(end)]);
    ylabel('\sigma [m^2]');
    xlabel('Wavelength [nm]')
    title('Cross Sections for Different Atmospheric Species')
    legend('O', 'N2', 'O2');
    xlim([0 40]);
    grid on;
        
end

%% READ IN LYRA Zr Channel Example Data

data = restore_idl([IDLfile_dir, LYRA_Zr_file]);
r_t = data.R_T;
n_des_col = data.N_COL_ARRAY*1E4;
extinction = data.RATIO_MEAS;

if plot_LYRA_Zr   
    LYRA_data = figure;
    plot(extinction, r_t-6371);
    xlabel('LYRA Extinction Ratio')
    ylabel('Altitude [km]')
    set(gca, 'XDir','reverse');     % reverse x values from high to low
end

%% READ IN the Response Function Data
%--------------------------------------------------------------------------
R_itp = zeros(num_filters, num_lambda);         % interpolated resp. funcs. for each filter at every lambda value from irradiance data

for i=1:num_filters-1
    data = restore_idl([IDLfile_dir, resp_func_file{i}]);
    
    if i == 1
        % flip the data so they are going from small - larger wavelengths
        lambda_filter = flip(data.W) ./10;          % wavelength of response functions [nm]
        R_Zr = flip(data.T_ZR);                     % responsitvity of Zirconium filter [electrons/photon]
        R_AlMg = flip(data.T_AL_MG);                % responsitvity of Aluminum-Magnesium filter [electrons/photon]

        % interpolate the response function to find R_AlMg and R_Zr at every lambda
        % value from the irradiance data
        R_itp(1,:) = interp1(lambda_filter, R_Zr, lambda);
        R_itp(2,:) = interp1(lambda_filter, R_AlMg, lambda);


        
    elseif i == 2
        lambda_filter_EUVM = data.RESP_WV;           % wavelength of response func [nm] this is the same for EUVMB
        R_EUVMA = data.RESP;                         % response function for EUVM A filter  [A/W]        
        R_itp(3,:) = interp1(lambda_filter_EUVM, R_EUVMA, lambda);
        
    elseif i == 3
        R_EUVMB = data.RESP;
        R_itp(4,:) = interp1(lambda_filter_EUVM, R_EUVMB, lambda);
    end
    
        % because the filter wavelengths are cut shorter than the irradiance
        % wavelengths, there are NaN values in the linear interpolated values. Get
        % rid of these by making them zero
        R_itp(isnan(R_itp)) = 0;
    
end

if plot_respFunc
    response_func = figure;
    
    for i=1:num_filters
        semilogy(lambda, R_itp(i,:), 'LineWidth', 2);
        hold on;
    end
    
    xlabel('Wavelength [nm]');
    xlim([lambda_filter(1) lambda_filter(end)]);
    ylabel('Responsivity [electrons/photon]');
    title('Responsivity of Foil Filters');
    legend('EUVOP Zr', 'EUVOP Al/Mg', 'EUVM Al/Nb/C', 'EUVM Al/Ti/C');
    xlim([0 40]);
    grid on;  
end

%% Convert Irradiance Into Photocurrent  
%--------------------------------------------------------------------------
diam = 5.5E-3;                  % detector diameter [m]
Ap = pi*(diam/2)^2;             % detector effective aperature [m^2]

E_photon = h*c./(lambda*1E-9);  % the energy of each photon (with wavelength lambda)
                                % from the irradiance data        [J/photon]
                          
photon_flux = zeros(size(irradiance));              

% find the photon flux for the irradiance data
num_days = length(dates);

for i=1:num_days
    for j=1:length(lambda)
        photon_flux(i,j) = irradiance(i,j) * Ap / E_photon(j);   % photon flux hitting the detector [photons/(s nm)]        
    end
end

% integrate over each response function separatly to find the total
% photocurrent 

% integrated photocurrents for each of the channels corresponding to each day
% of the irradiance data                                                    (numdays x 4)
I_photo = zeros(num_days, num_filters);                       % photocurrents [A]


for i=1:num_filters 
    % EUV OP resp funcs are in elec/photon
    if i < 3 
        I_photo(:,i) = trapz(lambda, (e .* photon_flux .*  R_itp(i,:)), 2);     % [A]
    % EUVM resp funcs are in A/W
    else
        I_photo(:,i) = trapz(lambda, (Ap .* irradiance .*  R_itp(i,:)), 2);           % [A]   
    end
             
end


% plot the photocurrents for each channel for year 2010 
if plot_photocurrent
    photo_curr = figure;
    for i=1:num_filters
        plot(I_photo((year == plotyear), i)./1E-9);
        hold on;
    end
    xlim([1 365]);
    xlabel('DOY');
    ylabel('Photocurrent [nA]');
    title(['Predicted Photocurrent at 1AU for ', num2str(plotyear)]);
    legend('EUVOP Zr', 'EUVOP Al/Mg', 'EUVM Al/Nb/C', 'EUVM Al/Ti/C');
    grid on;
end

%% Find LOS column number densites for every species i from MSIS
%--------------------------------------------------------------------------
 
N_LOS = zeros(num_z_occ, 3, num_dates);      % rows = altitude values, columns = [O, N2, O2] respectively, and page = dates

% for each occultation tangent altitude (row), consituent (column), and date (page),
% find the LOS column number dens
for k=1:num_dates
   fprintf('Finding LOS Column Number Densities for');
   DOY = doy(k,1)
   for j=1:3      
       if j == 1
           rho = rho_O;
           Constituent = 'O';
       elseif j == 2
           rho = rho_N2;
           Constituent = 'N2';
       else
           rho = rho_O2;
           Constituent = 'O2';
       end
       
       % run the tangent point heights from minz to second to last
       % point of occulation altitudes
       for i=1:num_z_occ
           z_h = z_occ(i);
           N_LOS(i, j, k) = N_integral(z_atm, rho(:,k), z_h, maxz_atm, z_sc, R);   % LOS column 
       end
   end    
end

% check the Column Number densities
if plot_col_densities
    col_dens = figure;
    plot(N_LOS(:, 1, which_date), z_occ/1000, N_LOS(:, 2, which_date), z_occ/1000, N_LOS(:, 3, which_date), z_occ/1000, 'linewidth', 2);
    legend('O', 'N2', 'O2', 'location', 'best');
    grid on;
    ylim([150 500])
    xlabel('LOS Integrated Column Number Density [1/m^2]');
    ylabel('Tangent Ray Point Alitude [km]');
    title('Integrated LOS Column Number Densities for 3/20/2010');
end

%% Calculate the Decay at every z_h for every lambda value using cross sections
%--------------------------------------------------------------------------
decay = zeros(num_z_occ, num_lambda, num_dates);    % exponential decay factor for every [altitude, lambda, and date]

for k=1:num_dates
    for i=1:num_z_occ
       for j=1:num_lambda
           
          val_O = sigma_O_itp(j) * N_LOS(i, 1, k);
          val_N2 = sigma_N2_itp(j) * N_LOS(i, 2, k);
          val_O2 = sigma_O2_itp(j) * N_LOS(i, 3, k);
          sum = val_O + val_N2 + val_O2;
          
          decay(i,j,k) = exp(-sum);         % value of exponential decay factor for every [altitude, lambda, and date]
          
       end
    end
end

% figure;
% cc = jet(num_lambda);
% for i=1:num_lambda
%     semilogx(decay(:, i, which_date), z_occ/1000, 'color', cc(i,:));
%     hold on;
% end
% xlim([1E-45 10]);
% xlabel('Decay Factor');
% ylabel('Tangent Ray Point Altitude [km]');

%% Calculate the Extinction Ratio
%--------------------------------------------------------------------------

% first calculate photocurrents at tangent ray height z_h when LOS is
% piercing atmosphere and signal is decaying
I_photo_z = zeros(num_z_occ, 2, num_filters); % rows = tangent ray height, col = date, page = filter  [A]   
S = zeros(num_z_occ, 2, num_filters);         % extinction ratio corresponding each tangent ray height [row] for each day [column] for every column [page]


% find the photon flux on each of the dates
date_index = zeros(num_dates,1);

days = [2010069, 2010265];

for i=1:num_dates
    date_index(i) = find(dates == days(i));
end

fprintf('Calculating Extinction Ratios...\n\nDONE! \n');

for k =1:num_filters
    for j = 1:num_dates
        for i=1:num_z_occ

            % Find the photocurrent at every altitude point after some extinction 
            if k < 3
                I_photo_z(i,j,k) = trapz(lambda,(e .* photon_flux(date_index(j),:) .* R_itp(k,:) .* decay(i, :, j) ), 2);
            %I_photo_AlMg_z(i,j) = trapz(lambda,(e .* photon_flux(date_index(j),:) .* (R_AlMg_itp)' .* decay(i, :, j) ), 2);

            else
                I_photo_z(i,j,k) = trapz(lambda,(Ap .* irradiance(date_index(j),:) .* R_itp(k,:) .* decay(i, :, j) ), 2);
                
            end
            
            % Extinction Ratio 
            S(i,j,k) = I_photo_z(i,j,k) / I_photo(date_index(j),k);             % rows = tangent point altitude, columns = dates, page = filter
            %S_AlMg(i,j) = I_photo_AlMg_z(i,j) / I_photo_AlMg(date_index(j));
        end
    end  
end

if plot_extinction
    ext_ratio = figure;
    
    for k=1:num_filters
        plot(S(:, which_date, k), z_occ/1000)
        hold on;
    end
    
    set(gca, 'XDir','reverse');     % reverse x values from high to low
    legend('EUVOP Zr', 'EUVOP Al/Mg', 'EUVM Al/Nb/C', 'EUVM Al/Ti/C');
    xlim([0 1]);
    xlabel('Normalized Extinction Ratio');
    ylabel('Tangent Ray Point Alitude [km]');
    title('Extinction Ratios 3/20/2010');
    grid on;
end   
%

%% ------------------- INVERSION THEORY -----------------------------------
%--------------------------------------------------------------------------
%--------------------- O + N2 Zirconium Filter ----------------------------

%% Calculate LOS column number density from "measured" LYRA Extinction Ratio 
% ***** N + O2 --> Zr filter only
%--------------------------------------------------------------------------
filter = 1;                     % 1 = EUVOP Zr, 2 = EUVOP Al/Mg, 3 = EUVM Al/Nb, 4 = EUVM Al/Ti/C
which_date = 1;                 % 1 = spring epnx, 2 = fall equinox
LYRA_ext = extinction;          % LYRA extinction ratios for every z_h below    
z_t_LYRA = r_t - R/1000;        % tangent ray height altitudes [km] for LYRA data
num_z_LYRA = length(z_t_LYRA);  % number of tangent ray points for LYRA data 
comb_CS = sigma_N2_itp + sigma_O_itp;   % combined cross section for N2 and O (for Zr channel) [m^2]
N_LOS_filt = zeros(num_z_LYRA,1);  % converged LOS colummn number dens for every z_h [1/m^2]

N_guess = 1E19;                  % initial guess for LOS column number density (highest altitude)

% LYRA data goes from low to high altitude, so this goes backwards (from high to low altitude)
for z = (num_z_LYRA:-1:1)
    S_measured = LYRA_ext(z);
              
    % N = LOS column number density to solve for
    f = @(N_LOS) func(N_LOS, S_measured, lambda, photon_flux(date_index(which_date),:), R_itp(filter,:), comb_CS);
    
    options = optimset('Display','notify');     % notify if a value doesnt converge
    % THE FACTOR OF TWO DOESN't MAKE SENSE
    N_LOS_filt(z) = 2*fzero(f, N_guess, options);
    
    N_guess = N_LOS_filt(z);               % use the next guess as the old one
end

% plot the inverted N_LOS data and compare with LYRA 
if plot_N_convergence_LYRA
    LOS_conv = figure;
    semilogx(N_LOS_filt, z_t_LYRA, 'r', n_des_col, z_t_LYRA, 'k')
    ylim([z_t_LYRA(1) z_t_LYRA(end)])
    xlabel('LOS Column Number Density [m^{-2}]');
    ylabel('Tangent Ray Point Altitude [km]');
    legend('N_{converged}', 'N_{LYRA}');
    title('Convered MSIS LOS Column Density VS LYRA Measurements for O + N2')
    grid on;
end
% 

%% Compute number density from LYRA measurements and compare with MSIS
% ***** N + O2 --> Zr filter only
%--------------------------------------------------------------------------

%----- create your L matrix -------
% rows => LOS tangent radius [m]
% cols => radial coordinate for number density, r_j 
 
% Since we are assuming a spherically symmetric atmosphere, the radial
% coordinate for number density is the radius of the atmospheric "shells"
% The elements of L_ij are zero for r_j values lower than the tangent ray
% height. This is because these lower shells don't attentuate the LOS ray.
% (it doesn't pass through them...)

r_MSIS = r_t*1000;          % LYRA tangent ray coordinates [m] from low to high
L = zeros(num_z_LYRA-1);      % LOS matrix [m]

for i = 1:num_z_LYRA-1            % LOS tangent     
   for j = i:num_z_LYRA-1       % radial coord       
       L(i,j) = sqrt( r_MSIS(j+1)^2 - r_MSIS(i)^2 ) - sqrt( r_MSIS(j)^2 - r_MSIS(i)^2 );
   end 
end

n_meas = (L' * L)^(-1) * L' * N_LOS_filt(1:end-1)/2;         % combined number denisty of O and N2 [1/m^3]
n_MSIS = rho_O(:,1) + rho_N2(:,1);                      % combined number density of O and N2 from MSIS [1/m^3]

if plot_inverted_dens_LYRA
    inversion_n = figure;
    semilogx(n_meas, z_t_LYRA(1:end-1))
    hold on;
    semilogx(n_MSIS, z_atm/1000);
    legend('LYRA Data', 'NRLMSIS-00')
    ylim([z_t_LYRA(1) z_t_LYRA(end-1)]);
    xlabel('Number Density [m^{-3}]');
    ylabel('Altitude [km]');
    title(' Combined Number Density of O and N2')
    grid on;
end

%% Compute number density from extinction ratio (simulated from MSIS) and compare with MSIS
% ***** N + O2 --> Zr filter only
%--------------------------------------------------------------------------

% Column Density First
% 1 = EUVOP Zr , 2 = EUVOP Al/Mg, 3 = EUVM Al/Nb, 4 = EUVM Al/Ti/C
% Zr = O and N2, Al/Mg = O2   

MSIS_z_t = z_occ/1000;            % tangent ray altitude [km]

comb_CS = sigma_N2_itp + sigma_O_itp;   % cross section for N2 and O for every lambda [m^2]
net_sigmas = [comb_CS, sigma_O_itp]; 
N_LOS_filt = zeros(num_z_occ,2);    % converged LOS colummn number dens for every z_h and filter [1/m^2]
N_guess = 1E19;                  % initial guess for LOS column number density (highest altitude)
z_cutoff = 100;                  % dont want to go below index of ~150 km

% MSIS z goes from low to high altitude, so this goes backwards (from high to low altitude)

for filt = 1:1      % 1 =  Zr (N2 + O)
    for z = (num_z_occ:-1:z_cutoff)
        S_measured = S(z, which_date, filt);        % extinction ratios computed from MSIS 

        % N_LOS_var = column number density to solve for
        f = @(N_LOS_var) func(N_LOS_var, S_measured, lambda, photon_flux(date_index(which_date),:), R_itp(filt,:), net_sigmas(:, filt));

        options = optimset('Display','notify');     % notify if a value doesnt converge
        % not sure why factor of 2 is there....
        N_LOS_filt(z, filt) = 2*fzero(f, N_guess, options);

        N_guess = N_LOS_filt(z,filt);                       % use the previous guess as the next one
    end
end

MSIS_N_LOS_comb = N_LOS(:, 1, which_date) + N_LOS(:, 2, which_date);  % combined column num dens. for O + N2

if plot_N_convergence_Zr
    figure;
    subplot(2,2,1);
    semilogx(MSIS_N_LOS_comb, MSIS_z_t, 'k', 'linewidth', 5);
    hold on;
    semilogx(N_LOS_filt(:,1), MSIS_z_t, '--g', 'linewidth', 2.5)
    ylim([MSIS_z_t(z_cutoff) MSIS_z_t(end)])
    xlabel('LOS Column Number Density [m^{-2}]');
    ylabel('Tangent Ray Point Altitude [km]');
    legend('N_{O+N2, MSIS}', 'N_{O+N2, Inverted}');
    title('O+N2 LOS Column Number Densities from Direct and Inverted MSIS Measurements')
    grid on;
end

% now number density
r_MSIS = z_occ + R;          % MSIS tangent ray coordinates [m] from low to high
L = zeros(num_z_occ-1);      % LOS matrix [m]


for i = 1:num_z_occ-1          % LOS tangent     
   for j = i:num_z_occ-1       % radial coord       
       L(i,j) = sqrt( r_MSIS(j+1)^2 - r_MSIS(i)^2 ) - sqrt( r_MSIS(j)^2 - r_MSIS(i)^2 );
   end 
end

n_meas = (L' * L)^(-1) * L' * N_LOS_filt((1:end-1), 1)/1.7;    % combined number denisty of O and N2 [1/m^3]
n_MSIS = rho_O(:,1) + rho_N2(:,1);                           % combined number density of O and N2 from MSIS [1/m^3]

if plot_inverted_dens_Zr
    subplot(2,2,2);
    semilogx(n_MSIS, z_atm/1000, 'linewidth', 2);
    hold on;
    semilogx(n_meas, MSIS_z_t(1:end-1), '-.', 'linewidth', 2);
    legend('n_{O+N2, MSIS}', 'n_{O+N2, Inverted}')
    ylim([MSIS_z_t(z_cutoff) MSIS_z_t(end-1)]);
    xlabel('Number Density [m^{-3}]');
    ylabel('Altitude [km]');
    title('O+N2 Number Densities from Direct and Inverted MSIS Measurements')
    grid on;
end

% --------------------- O/N2 Zr vs. Al/Mg Filter ----------------------------
x0 = [2E19, 3E17];                  % initial LOS number dens guesses for TOA 1 = Zr, 2 = Al/Mg 
N_LOS_filt2 = zeros(num_z_occ,2);
options = optimoptions('fsolve', 'FiniteDifferenceType', 'central', 'OptimalityTolerance', 1E-8, 'Display', 'off');

fprintf('\nCalculating LOS Number Densities for O/N2 Ratio... \n\n');

for z = (num_z_occ:-1:z_cutoff) 
    params1 = [S(z, which_date, 1), S(z, which_date, 2)];    
    % make everything in params2 a column vector...
    params2 = [lambda, photon_flux(date_index(which_date),:)', R_itp(1,:)', R_itp(2,:)', sigma_O_itp, sigma_N2_itp];
        
    % the variable x will be the N_LOS values for the different filters
    % x(1) = Zr, x(2) = Al/Mg
    
    % in order to get fsolve to work correctly, we must scale x
    % accordingly.
    scaleFactor = 1E20;   
    x0 = x0 ./ scaleFactor;
    
    f = @(x)func2(x, params1, params2, scaleFactor);
    [xsol, fval] = fsolve(f, x0, options);
    
    N_LOS_filt2(z,:) = xsol .* scaleFactor;
    
    x0 = N_LOS_filt2(z,:);
end

% get rid of outliers and replace 
N_LOS_filt2 = filloutliers(N_LOS_filt2, 'linear', 'movmedian', 3, 'ThresholdFactor', 2 );

fprintf('DONE! \n\n');

if plot_N_convergence_AlMg
    subplot(2,2,3);
    % actual MSIS output
    semilogx(N_LOS(:, 1, which_date), MSIS_z_t, 'color', [0.4940, 0.1840, 0.5560], 'linewidth', 5)
    hold on;
    semilogx(N_LOS(:, 2, which_date), MSIS_z_t, 'k', 'linewidth', 5) 
    hold on;    
    % inverted values
    semilogx(N_LOS_filt2(:,1), MSIS_z_t, '--y', 'linewidth', 2.5)
    hold on;
    semilogx(N_LOS_filt2(:,2), MSIS_z_t, '--g', 'linewidth', 2.5);
    
    ylim([MSIS_z_t(z_cutoff) MSIS_z_t(end)])
    xlabel('LOS Column Number Density [m^{-2}]');
    ylabel('Tangent Ray Point Altitude [km]');
    legend(' N_{O, MSIS}', 'N_{N2, MSIS}', 'N_{O, inverted}', 'N_{N2, inverted}');
    title('O and N2 LOS Column Number Densities from Direct and Inverted MSIS Measurements')
    grid on;
end


% now compute number density as a function of z
% we already have the L matrix and r_MSIS (the tangent ray coords) from the
% Zr filter, so we can use these for the AlMg filter as well.

n_meas2 = (L' * L)^(-1) * L' * N_LOS_filt2((1:end-1), :)/1.7; % number densities of O and N2 [1/m^3]
n_MSIS2 = [rho_O(:,1)  rho_N2(:,1)];                      % number densities of O and N2 from MSIS [1/m^3]

% remove and replace outliers
n_meas2 = filloutliers(n_meas2, 'linear', 'movmedian', 5, 'ThresholdFactor', 2 );

if plot_inverted_dens_AlMg
    subplot(2,2,4);
    % actual MSIS output
    semilogx(n_MSIS2(:, 1), z_atm/1000, 'linewidth', 2)
    hold on;
    semilogx(n_MSIS2(:, 2), z_atm/1000, 'linewidth', 2) 
    hold on;    
    % inverted values
    semilogx(n_meas2(:,1), MSIS_z_t(1:end-1), '-.', 'linewidth', 2)
    hold on;
    semilogx(n_meas2(:,2), MSIS_z_t(1:end-1), '-.', 'linewidth', 2);

    legend(' n_{O, MSIS}', 'n_{N2, MSIS}', 'n_{O, inverted}', 'n_{N2, inverted}');
    ylim([MSIS_z_t(z_cutoff) MSIS_z_t(end-1)]);
    xlabel('Number Density [m^{-3}]');
    ylabel('Altitude [km]');
    title('O and N2 Number Densities from Direct and Inverted MSIS Measurements')
    grid on;
end



































