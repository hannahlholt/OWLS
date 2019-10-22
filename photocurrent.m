
close all;


% EDIT THESE.....
%-------------------------------------------------------------------------
IDLfile_dir = '/Users/hannahholt/Documents/IDLfiles/';  %directory with .sav files
irradiance_file = 'FISM_EVE_daily_merged_earth_v8r1_01_2007_on.sav';
resp_func_file = 'EUV_OP_resp_function.sav';
cross_sec_file = 'photon_cross_sections.sav';
LYRA_Zr_file = 'N_col_o_2014-12-20_41.9064_km.sav';

plot_irradiance = false;            % To plot an example irradiance
plot_respFunc = false;              % To plot the response functions for each channels
plot_cross_section = false;         % To plot the photon cross sections
plot_LYRA_Zr = false;                % To plot the example data set from the LYRA Zr channel
plot_photocurrent = false;           % To plot the photocurrents
plot_densities = false;             % To plot MSIS mass density outputs
plot_col_densities = false;         % To plot the integrate LOS column number densities
plot_extinction = false;             % To plot the extinction ratio for the equinoxes

plotyear = 2010;                    % which year do you want to plot irradiances
plotlambda = 10;                    % which wavelength do you want to plot [nm]
which_date = 1;                     % which date you want to plot (1 = spring eqnx, 2 = fall eqnx)

R_sun = 6.957E8;                    % Solar radius [m]
b = 1.49597870700E11;               % avg distance from Earth to Sun (m)
zeta_max = atan(R_sun/b);           % Ang from center of Sun to edge of disk (AS SEEN FROM EARTH!) [rad]
N_zeta_pts = 100;
dzeta = zeta_max/N_zeta_pts;        
%-------------------------------------------------------------------------


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

%% READ IN the Response Function Data
%--------------------------------------------------------------------------

data = restore_idl([IDLfile_dir, resp_func_file]);

% flip the data so they are going from small - larger wavelengths
lambda_filter = flip(data.W) ./10;          % wavelength of response functions [nm]
R_Zr = flip(data.T_ZR);                     % responsitvity of Zirconium filter [electrons/photon]
R_AlMg = flip(data.T_AL_MG);                % responsitvity of Aluminum-Magnesium filter [electrons/photon]

% interpolate the response function to find R_AlMg and R_Zr at every lambda
% value from the irradiance data
R_Zr_itp = interp1(lambda_filter, R_Zr, lambda);
R_AlMg_itp = interp1(lambda_filter, R_AlMg, lambda);

% because the filter wavelengths are cut shorter than the irradiance
% wavelengths, there are NaN values in the linear interpolated values. Get
% rid of these by making them zero
R_Zr_itp(isnan(R_Zr_itp)) = 0;
R_AlMg_itp(isnan(R_AlMg_itp)) = 0;

if plot_respFunc
    response_func = figure;
    semilogy(lambda_filter, R_Zr, lambda_filter, R_AlMg, 'LineWidth', 2);
    %semilogy(lambda, R_Zr_itp, lambda, R_AlMg_itp, 'LineWidth', 2);
    xlabel('Wavelength [nm]');
    xlim([lambda_filter(1) lambda_filter(end)]);
    ylabel('Responsivity [electrons/photon]');
    title('Responsivity of EUV-OP Foil Filters');
    legend('Zr', 'Al-Mg')
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
    subplot(3,1,1);
    plot(lambda_CS, sigma_O);
    xlim([lambda_CS(1), lambda_CS(end)]);
    ylabel('\sigma [m^2]');
    title('Cross Sections for Different Atmospheric Species')
    legend('O');
    grid on;
    
    subplot(3,1,2)
    plot(lambda_CS, sigma_N2, 'r');
    xlim([lambda_CS(1), lambda_CS(end)]);
    ylabel('\sigma [m^2]');
    legend('N2');
    grid on;
    
    subplot(3,1,3)
    plot(lambda_CS, sigma_O2, 'g');
    xlim([lambda_CS(1), lambda_CS(end)]);
    xlabel('Wavelength [nm]');
    ylabel('\sigma [m^2]');
    legend('O2');
    grid on;
end

%% READ IN LYRA Zr Channel Example Data

data = restore_idl([IDLfile_dir, LYRA_Zr_file]);

r_t = data.R_T;
n_des_col = data.N_COL_ARRAY*1E4;
extinction = data.RATIO_MEAS;
%plot(extinction, r_t-6371)
%xlabel('column dens')
%ylabel('Altitude [m]')


%% Convert Irradiance Into Photocurrent  
%--------------------------------------------------------------------------
diam = 5.5E-3;                  % detector diameter [m]
Ap = pi*(diam/2)^2;             % detector effective aperature [m^2]

E_photon = h*c./(lambda*1E-9);  % the energy of each photon (with wavelength lambda)
                                % from the irradiance data [J/photon]
                          
photon_flux = zeros(size(irradiance));              

% find the photon flux for the irradiance data
for i=1:length(dates)
    for j=1:length(lambda)
        photon_flux(i,j) = irradiance(i,j) * Ap / E_photon(j);   % photon flux hitting the detector [photons/(s nm)]        
    end
end

% integrate over each response function separatly to find the total
% photocurrent 

% integrated photocurrents for each of the channels corresponding to each day
% of the irradiance data                                                    (nx1)
I_photo_Zr = trapz(lambda,(e .* photon_flux .* (R_Zr_itp)'),2);         % [A]
I_photo_AlMg = trapz(lambda,(e .* photon_flux .* (R_AlMg_itp)'),2);     % [A]

% plot the photocurrents for each channel for year 2010 
if plot_photocurrent
    photo_curr = figure;
%     plot(plotdoy, irradiance((year == plotyear), index));
    plot(plotdoy, I_photo_Zr(year == plotyear)./1E-9);
    hold on;
    plot(plotdoy, I_photo_AlMg(year == plotyear)./1E-9);
    xlim([1 365]);
    xlabel('DOY');
    ylabel('Photocurrent [nA]');
    title(['Predicted Photocurrent at 1AU for ', num2str(plotyear)]);
    legend('Zr', 'AlMg');
    grid on;
end

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
    ylabel('Altitude [km]');
    title('NRLMSIS-00 Mass Densities for 3/20/2010');
end

%% Find LOS column number densites for every species i
%--------------------------------------------------------------------------
 
N_LOS = zeros(num_z_occ, 3, num_dates);      % rows = altitude values, columns = [O, N2, O2] respectively, and page = dates

% for each occultation tangent altitude (row), consituent (column), and date (page),
% find the LOS column number dens
for k=1:num_dates
    fprintf('Finding LOS Column Number Densities for');
   DOY = doy(1,k)
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
I_photo_Zr_z = zeros(num_z_occ, 2);         % rows = tangent ray height, col = date  [A]      
I_photo_AlMg_z = zeros(num_z_occ, 2);
S_Zr = zeros(num_z_occ, 2);                 % extinction ratio corresponding each tangent ray height [row] for each day [column]
S_AlMg = zeros(num_z_occ, 2);


% find the photon flux on each of the dates
date_index = zeros(num_dates,1);

days = [2010069, 2010265];

for i=1:num_dates
    date_index(i) = find(dates == days(i));
end

fprintf('Calculating Extinction Ratios...\n\nDONE! \n');

for k = 1:num_dates
    for i=1:num_z_occ
        
        % Find the photocurrent at every altitude point after some extinction 
        I_photo_Zr_z(i,k) = trapz(lambda,(e .* photon_flux(date_index(k),:) .* (R_Zr_itp)' .* decay(i, :, k) ), 2);
        I_photo_AlMg_z(i,k) = trapz(lambda,(e .* photon_flux(date_index(k),:) .* (R_AlMg_itp)' .* decay(i, :, k) ), 2);
        
        % Extinction Ratio 
        S_Zr(i,k) = I_photo_Zr_z(i,k) / I_photo_Zr(date_index(k));             % rows = tangent point altitude, columns = dates
        S_AlMg(i,k) = I_photo_AlMg_z(i,k) / I_photo_AlMg(date_index(k));
    end
    
end

if plot_extinction
    ext_ratio = figure;
    subplot(1,2,1)
    plot(S_Zr(:, which_date), z_occ/1000, S_AlMg(:, which_date), z_occ/1000)
    set(gca, 'XDir','reverse');     % reverse x values from high to low
    legend('Zr', 'AlMg')
    xlim([0 1]);
    xlabel('Normalized Extinction Ratio');
    ylabel('Tangent Ray Point Alitude [km]');
    title('Extinction Ratios 3/20/2010');
    grid on;
    
    
    subplot(1,2,2)
    plot(S_Zr(:, which_date+1), z_occ/1000, S_AlMg(:, which_date+1), z_occ/1000)
    set(gca, 'XDir','reverse')
    legend('Zr', 'AlMg');
    xlim([0 1]);
    xlabel('Normalized Extinction Ratio');
    ylabel('Tangent Ray Point Alitude [km]');
    title('Extinction Ratios for 9/22/2010');
    grid on;
    
end

%




