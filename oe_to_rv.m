% 
% Calculates position and velocity vectors from orbital elements.
% INPUTS
%  a_     - semimajor axis                        (km)
%  e_     - eccentricity                          (unitless)
%  i_     - inclination                           (degrees)
%  omega_ - right ascension of the ascending node (degrees)
%  w_     - argument of periapse                  (degrees)
%  M_     - mean anomaly                          (degrees)
%  mu_    - gravitational parameter               (km^3/s^2)
% OUTPUTS
%  r_ - position vector (km)
%  v_ - velocity vector (km/s)
function[r_, v_] = oe_to_rv(a_, e_, i_, omega_, w_, M_, mu_)

    %convert values to radians
    i_ = i_*pi/180;
    omega_ = omega_*pi/180;
    w_ = w_*pi/180;
    M_ = M_*pi/180;

    %calculate eccentric anomaly
    syms E_;
    Keplerseqn = M_ == E_ - e_*sin(E_);
    E_ = vpasolve(Keplerseqn, E_);

    %calculate true anomaly
    f_ = vpa(2*atan(sqrt((1+e_)/(1-e_))*tan(E_/2)));

    %find angular momentum
    h_ = sqrt(mu_*a_*(1-e_^2));

    %find theta
    theta = vpa(w_ + f_);

    %calculate r
    r = a_*(1 - e_^2)/(1 + e_*cos(f_));

    %find r vector
    ri = cos(theta)*cos(omega_) - cos(i_)*sin(omega_)*sin(theta);
    rj = cos(theta)*sin(omega_) + cos(i_)*cos(omega_)*sin(theta);
    rk = sin(i_)*sin(theta);
    r_ = vpa(r*[ri, rj, rk]); %km

    %find v vector
    vi = -(mu_/h_)*(cos(omega_)*(sin(theta) + e_*sin(w_)) + sin(omega_)*(cos(theta) + e_*cos(w_))*cos(i_));
    vj = -(mu_/h_)*(sin(omega_)*(sin(theta) + e_*sin(w_)) - cos(omega_)*(cos(theta) + e_*cos(w_))*cos(i_));
    vk = (mu_/h_)*((cos(theta) + e_*cos(w_))*sin(i_));
    v_ = vpa([vi, vj, vk]); %km/s

    r_ = double(r_);
    v_ = double(v_);

end