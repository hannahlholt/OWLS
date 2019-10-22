function [N_LOS] = N_integral(z, n_i, R)
% N_integral is a function that computes the LOS column number density of a
% species i when doing an occulation measurement. Returns a double that
% corresponds to the LOS column number density at for a given orbital
% altitude z_sc and altitude of the vector from Earth's surface to the LOS
% tangent vector. Assumes z and n_i have same altitude spacing and length

% ------- INPUTS -----------
% z         - vector of altitude values for atmosphere      [m]        (nx1)                       
% n_i       - vector of number densities for species i      [1/m^3]    (nx1)  
%             for all altitude values of z              
% R         - radius of Earth (volumetric)                  [m]        (single)

% ------- RETURNS -----------
% N_LOS     - LOS column number density for species i       [m^2]      (single)
%             for given z_h and z_sc               

r = z + R;          % radial value of tangent point;
num_r = length(r);

L = zeros(num_r-1);

for i = 1:num_r-1              % LOS tangent     
   for j = i:num_r-1           % radial coord       
       L(i,j) = sqrt( r(j+1)^2 - r(i)^2 ) - sqrt( r(j)^2 - r(i)^2 );
   end
end

% multiply by n to get N_LOS at every z;
N_LOS = 2 * n_i(1:end-1)' * L';

end

