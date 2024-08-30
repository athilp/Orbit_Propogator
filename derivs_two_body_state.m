function [f] = derivs_two_body_state(state, params, v)
c_D = params(1);
a = params(2);
m = params(3);
J2 = params(4);
J3 = params(5);
mu = params(6);
R_E = params(7);
td = params(8);

    %Initialize the derivative vector
    
    f = zeros(6,1);

    %States
    x = state(1);
    y = state(2);
    z = state(3);
    xd = state(4);
    yd = state(5);
    zd = state(6);
    
    r = sqrt(x^2+y^2+z^2);

    %Get the density at altitude
    h = r-R_E;
    [rho0_kgm3, h0, H] = getDensityParams(h); 

    %Convert rho0 from kg/m^3 to kg/km^3
    rho0_kgkm3 = rho0_kgm3*10^9;
   

    %Equations of Motion
    f(1) = xd;

    f(2) = yd;

    f(3) = zd;

    f(4) = (3*J2*R_E^2*mu*x*z^2)/(x^2 + y^2 + z^2)^(7/2) - (J3*R_E^3*mu*((3*x*z)/(x^2 + y^2 + z^2)^(3/2) - (15*x*z^3)/(x^2 + y^2 + z^2)^(5/2)))/(2*(x^2 + y^2 + z^2)^2) - (mu*x)/(x^2 + y^2 + z^2)^(3/2) + (9*J2*R_E^2*mu*x*(z^2/(x^2 + y^2 + z^2) - 1/3))/(2*(x^2 + y^2 + z^2)^(5/2)) - (2*J3*R_E^3*mu*x*((3*z)/(x^2 + y^2 + z^2)^(1/2) - (5*z^3)/(x^2 + y^2 + z^2)^(3/2)))/(x^2 + y^2 + z^2)^3 + ((-c_D*a)/(2*m))*(rho0_kgkm3*exp((h0 + R_E - sqrt(x^2+y^2+z^2))/H))*sqrt((xd+td*y)^2 + (yd-td*x)^2 +zd^2)*(xd+td*y) +v(1);


    f(5) = (3*J2*R_E^2*mu*y*z^2)/(x^2 + y^2 + z^2)^(7/2) - (J3*R_E^3*mu*((3*y*z)/(x^2 + y^2 + z^2)^(3/2) - (15*y*z^3)/(x^2 + y^2 + z^2)^(5/2)))/(2*(x^2 + y^2 + z^2)^2) - (mu*y)/(x^2 + y^2 + z^2)^(3/2) + (9*J2*R_E^2*mu*y*(z^2/(x^2 + y^2 + z^2) - 1/3))/(2*(x^2 + y^2 + z^2)^(5/2)) - (2*J3*R_E^3*mu*y*((3*z)/(x^2 + y^2 + z^2)^(1/2) - (5*z^3)/(x^2 + y^2 + z^2)^(3/2)))/(x^2 + y^2 + z^2)^3 + ((-c_D*a)/(2*m))*(rho0_kgkm3*exp((h0 + R_E - sqrt(x^2+y^2+z^2))/H))*sqrt((xd+td*y)^2 + (yd-td*x)^2 +zd^2)*(yd-td*x) + v(2);


    f(6) = (J3*R_E^3*mu*(3/(x^2 + y^2 + z^2)^(1/2) - (18*z^2)/(x^2 + y^2 + z^2)^(3/2) + (15*z^4)/(x^2 + y^2 + z^2)^(5/2)))/(2*(x^2 + y^2 + z^2)^2) - (3*J2*R_E^2*mu*((2*z)/(x^2 + y^2 + z^2) - (2*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)^(3/2)) - (mu*z)/(x^2 + y^2 + z^2)^(3/2) + (9*J2*R_E^2*mu*z*(z^2/(x^2 + y^2 + z^2) - 1/3))/(2*(x^2 + y^2 + z^2)^(5/2)) - (2*J3*R_E^3*mu*z*((3*z)/(x^2 + y^2 + z^2)^(1/2) - (5*z^3)/(x^2 + y^2 + z^2)^(3/2)))/(x^2 + y^2 + z^2)^3 + ((-c_D*a)/(2*m))*(rho0_kgkm3*exp((h0 + R_E - sqrt(x^2+y^2+z^2))/H))*sqrt((xd+td*y)^2 + (yd-td*x)^2 +zd^2)*(zd)+v(3);


    