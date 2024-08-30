function [statedot] = derivs_two_body(full_state, params, v)
    

c_D = params(1);
a = params(2);
m = params(3);
J2 = params(4);
J3 = params(5);
mu = params(6);
R_E = params(7);
td = params(8);


    %Initialize the derivative vector
    
    statedot = zeros(42,1);

    %States
    x = full_state(1);
    y = full_state(2);
    z = full_state(3);
    xd = full_state(4);
    yd = full_state(5);
    zd = full_state(6);
    STM = full_state(7:42);
    STM = reshape(STM,6,6);

    state = full_state(1:6);

    %Useful expressions/definitions
    r = sqrt(x^2+y^2+z^2);

    %Get the density at altitude
    h = r-R_E;
    [rho0_kgm3, h0, H] = getDensityParams(h); 

    %Convert rho0 from kg/m^3 to kg/km^3
    rho0_kgkm3 = rho0_kgm3*10^9;
   

    %Equations of Motion
    statedot(1) = xd;

    statedot(2) = yd;

    statedot(3) = zd;

    statedot(4) = (3*J2*R_E^2*mu*x*z^2)/(x^2 + y^2 + z^2)^(7/2) - (J3*R_E^3*mu*((3*x*z)/(x^2 + y^2 + z^2)^(3/2) - (15*x*z^3)/(x^2 + y^2 + z^2)^(5/2)))/(2*(x^2 + y^2 + z^2)^2) - (mu*x)/(x^2 + y^2 + z^2)^(3/2) + (9*J2*R_E^2*mu*x*(z^2/(x^2 + y^2 + z^2) - 1/3))/(2*(x^2 + y^2 + z^2)^(5/2)) - (2*J3*R_E^3*mu*x*((3*z)/(x^2 + y^2 + z^2)^(1/2) - (5*z^3)/(x^2 + y^2 + z^2)^(3/2)))/(x^2 + y^2 + z^2)^3 + ((-c_D*a)/(2*m))*(rho0_kgkm3*exp((h0 + R_E - sqrt(x^2+y^2+z^2))/H))*sqrt((xd+td*y)^2 + (yd-td*x)^2 +zd^2)*(xd+td*y) + v(1);


    statedot(5) = (3*J2*R_E^2*mu*y*z^2)/(x^2 + y^2 + z^2)^(7/2) - (J3*R_E^3*mu*((3*y*z)/(x^2 + y^2 + z^2)^(3/2) - (15*y*z^3)/(x^2 + y^2 + z^2)^(5/2)))/(2*(x^2 + y^2 + z^2)^2) - (mu*y)/(x^2 + y^2 + z^2)^(3/2) + (9*J2*R_E^2*mu*y*(z^2/(x^2 + y^2 + z^2) - 1/3))/(2*(x^2 + y^2 + z^2)^(5/2)) - (2*J3*R_E^3*mu*y*((3*z)/(x^2 + y^2 + z^2)^(1/2) - (5*z^3)/(x^2 + y^2 + z^2)^(3/2)))/(x^2 + y^2 + z^2)^3 + ((-c_D*a)/(2*m))*(rho0_kgkm3*exp((h0 + R_E - sqrt(x^2+y^2+z^2))/H))*sqrt((xd+td*y)^2 + (yd-td*x)^2 +zd^2)*(yd-td*x) + v(2);


    statedot(6) = (J3*R_E^3*mu*(3/(x^2 + y^2 + z^2)^(1/2) - (18*z^2)/(x^2 + y^2 + z^2)^(3/2) + (15*z^4)/(x^2 + y^2 + z^2)^(5/2)))/(2*(x^2 + y^2 + z^2)^2) - (3*J2*R_E^2*mu*((2*z)/(x^2 + y^2 + z^2) - (2*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)^(3/2)) - (mu*z)/(x^2 + y^2 + z^2)^(3/2) + (9*J2*R_E^2*mu*z*(z^2/(x^2 + y^2 + z^2) - 1/3))/(2*(x^2 + y^2 + z^2)^(5/2)) - (2*J3*R_E^3*mu*z*((3*z)/(x^2 + y^2 + z^2)^(1/2) - (5*z^3)/(x^2 + y^2 + z^2)^(3/2)))/(x^2 + y^2 + z^2)^3 + ((-c_D*a)/(2*m))*(rho0_kgkm3*exp((h0 + R_E - sqrt(x^2+y^2+z^2))/H))*sqrt((xd+td*y)^2 + (yd-td*x)^2 +zd^2)*(zd)+ v(3);


    
    %Analytically derivied
%   statedot(4) = -x*mu/r^3+((J2*mu*R_E^2)/2)*(15*x*z^2/r^7-3*x/r^5)+((J3*mu*R_E^3)/2)*(35*x*z^3/r^9-15*x*z/r^7) -(c_D*a*pho_A/(2*m))*(xd+td*y)*v;
%   statedot(5) = -y*mu/r^3+((J2*mu*R_E^2)/2)*(15*y*z^2/r^7-3*y/r^5)+((J3*mu*R_E^3)/2)*(35*y*z^3/r^9-15*y*z/r^7) -(c_D*a*pho_A/(2*m))*(yd-td*x)*v;
%   statedot(6) = -z*mu/r^3-((J2*mu*R_E^2)/2)*(3*z^2/r^5+(6*x^2*z+6*y^2*z-9*z^3)/r^7)+((J3*mu*R_E^3)/2)*((3*x^2+3*y^2-12*z^2)/r^7-(15*x^2*z^2+15*y^2*z^2-20*z^4)/r^9) -(c_D*a*pho_A/(2*m))*zd*v;
    

    %Define the A matrix and evaluate at nonlinear trajectory
    
    %Integration of STM
    A = computeA(state, params);
    
    STM_dot = A*STM;

    %Reshape the vectors back to a nx1
    statedot_7_42 = reshape(STM_dot, 36,1);
    
    statedot = [statedot(1);statedot(2);statedot(3);statedot(4);statedot(5);statedot(6);statedot_7_42];
    
end