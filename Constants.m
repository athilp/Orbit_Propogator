

%Base case
c_D = 2;   %N/A
a = 0.0000036;   %m^2
m = 1350;  %kg
J2 = 0.0010826267;
J3 = -0.0000025327;
mu = 398600.4415;  %km^3/s^2
R_E = 6378.1363;   %km
thetadot = 7.29211585530066*10^-5;  %rad/s


%1 - Coefficient of  Drag, 2 - Area, 3 - Mass, 4 - J2, 5 - J3, 6 - mu, 
% 7 - R_E,, 8 - thetadot
params = [c_D; a; m;J2;J3;mu;R_E;thetadot];