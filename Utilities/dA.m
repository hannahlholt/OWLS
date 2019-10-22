function [ dA_value ] = dA( zeta_max, zeta_i, b, dzeta )
% Function to calculate the area element dA on the sun. zetas are in
% radians

dA_value = -2 * sqrt( tan(zeta_max)^2 - tan(zeta_i)^2 )* b^2 * sec(zeta_i)^2 * dzeta;


end

