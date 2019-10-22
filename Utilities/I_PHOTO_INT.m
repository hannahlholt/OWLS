function [ S ] = I_PHOTO_INT(N_LOS, lambda, photon_flux, R, sigma )
% I_PHOTO_INT returns the extinction ratio S at a given tangent ray point
% (using integration). S = decayed signal/TOA signal

e = 1.60217662e-19;                             % elementary charge [C]

% ------- INPUTS -----------
% lambda        - vector of wavelengths  to integrate over          [nm]    (nx1)    
% photon_flux   - vector of photon fluxes for a specific z_h for every 
%                 wavelength (NOT multiplied by response func)      [photons/(s nm)]  (nx1)                       
% R             - vector of response function vs. lambda            []      (nx1)  
% sigma         - vector of species cross-sec vs. lambda            [m^2]   (nx1)
% N_LOS         - LOS column number density for a specific species  [1/m^2] (single)

% ------- RETURNS -----------
% S             - extinction ratio for the specific z_h             []      (single)
%                            
 
% Top of atmosphere photocurrent
denom = trapz(lambda, (e .* photon_flux .*  R), 2);                   % [A]

% exponential decay factor for every value of lambda given by sigma (row vector)
decay_fac = (exp(-double(sigma) .* N_LOS))';  

% extinguished photocurrent 
num = trapz(lambda, (e .* photon_flux .*  R .* decay_fac), 2);        % [A]

S = double(num/denom);

end

