function [value] = func(N_LOS, S_measured, lambda, photon_flux, R, sigma)      
    S_calc = I_PHOTO_INT(N_LOS, lambda, photon_flux, R, sigma);
    value = S_calc - S_measured;
end

