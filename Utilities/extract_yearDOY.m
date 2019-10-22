function [ year, doy ] = extract_yearDOY( dates )
% This function takes in a vector of doubles corresponding to the year and doy 
% [yyyydoy] for consecutive times. Converts into two seperate vectors the length of
% the input, one with the years, and one with the doy.

year = num2str(dates); 
year = str2num(year(:,1:4));
doy = num2str(dates); 
doy = str2num(doy(:,5:end));  

end

