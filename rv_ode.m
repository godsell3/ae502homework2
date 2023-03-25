% 
% ODE used to solve for an orbit given a state vector of r and v vectors
% and gravitational parameter. Incorporates the J2 perturbation.
% INPUTS
%  state - condition vector containing position and velocity at a specific 
%          time [r1 r2 r3 v1 v2 v3] (km km km km/s km/s km/s)
%  mu    - gravitational parameter  (km^3/s^2)
%  J2    - J2 perturbation value    (unitless)
%  R     - radius of orbited body   (km)
% OUTPUTS
%  endcond - condition vector after taking a time step [r1 r2 r3 v1 v2 v3] 
%            (km km km km/s km/s km/s)
function[endcond] = rv_ode(~, state, mu, J2, R)
    rmag = norm(state(1:3));
    coeff = -mu/rmag^3;
    
    x = state(1);
    y = state(2);
    z = state(3);
    pcoeff = (3/2)*(J2*mu*R^2)/(rmag^4);
    pJ2(1) = pcoeff*(x/rmag)*(5*(z/rmag)^2 - 1);
    pJ2(2) = pcoeff*(y/rmag)*(5*(z/rmag)^2 - 1);
    pJ2(3) = pcoeff*(z/rmag)*(5*(z/rmag)^2 - 3);
    
    endcond = [state(4); state(5); state(6); coeff*state(1)+pJ2(1); coeff*state(2)+pJ2(2); coeff*state(3)+pJ2(3)];
end