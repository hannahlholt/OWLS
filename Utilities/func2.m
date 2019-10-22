function [eqns] = func2(x, params1, params2, scaleFactor)
% This function is used to calculate the number densities for the  Zr and Al/Mg
% filters, which decouples the N2 and O ratios and solves them
% simultaneously. 
% cross sections are kept separate. 

    e = 1.60217662e-19;                             % elementary charge [C]

    % Unscale x
    x = x .* scaleFactor;

    S_meas_Zr = params1(1);
    S_meas_AlMg = params1(2);

    lambda = params2(:,1);
    photon_flux = params2(:,2);
    R_Zr = params2(:,3);
    R_AlMg = params2(:,4);
    sigmaO = params2(:,5);
    sigmaN2 = params2(:,6);
    
   
    
    S_calc_Zr = I_PHOTO_INT2(x, lambda, photon_flux, R_Zr, sigmaO, sigmaN2);
    S_calc_AlMg = I_PHOTO_INT2(x, lambda, photon_flux, R_AlMg, sigmaO, sigmaN2);  
    
        
    % maybe switch these....    
    eqn1 = S_calc_Zr - S_meas_Zr; 
    eqn2 = S_calc_AlMg - S_meas_AlMg;
    eqns = [eqn1; eqn2];
   
end

