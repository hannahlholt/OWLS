function [ S ] = I_PHOTO_INT2(x, lambda, photon_flux, R, sigmaO, sigmaN2)
% I_PHOTO_INT returns the extinction ratio S at a given tangent ray point
% (using integration). S = decayed signal/TOA signal

e = 1.60217662e-19;                             % elementary charge [C]

% ------- INPUTS -----------
% lambda        - vector of wavelengths  to integrate over          [nm]    (nx1)    
% photon_flux   - vector of photon fluxes for a specific z_h for every 
%                 wavelength (NOT multiplied by response func)      [photons/(s nm)]  (nx1)                       
% R             - vector of response function vs. lambda            []      (nx1)  
% sigma_O       - vector of species cross-sec vs. lambda            [m^2]   (nx1)
% sigma_N2 
% N_LOS         - LOS column number density for a specific species  [1/m^2] (single)

% ------- RETURNS -----------
% S             - extinction ratio for the specific z_h             []      (single)
%                            
 

% Top of atmosphere photocurrent
denom = trapz(lambda, (e .* photon_flux .*  R), 1);                   % [A] dim = 1 to integrate the column. returns a row vector

% exponential decay factor for every value of lambda given by sigma (row vector)
decay_fac = exp(-( (sigmaO .* x(1)) + (sigmaN2 .* x(2) ) ) );  

% extinguished photocurrent 
num = trapz(lambda, (e .* photon_flux .*  R .* decay_fac), 1);        % [A]

S = double(num./denom);


end

