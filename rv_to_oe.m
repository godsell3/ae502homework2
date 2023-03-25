% 
% Calculates orbital elements from position and velocity vectors.
% INPUTS
%  r_  - position vector         (km)
%  v_  - velocity vector         (km/s)
%  mu_ - gravitational parameter (km^3/s^2)
% OUTPUTS
%  a_    - semimajor axis                        (km)
%  e_    - eccentricity                          (unitless)
%  i_    - inclination                           (degrees)
%  RAAN_ - right ascension of the ascending node (degrees)
%  w_    - argument of periapse                  (degrees)
%  f_    - true anomaly                          (degrees)
function[a_, e_, i_, RAAN_, w_, f_] = rv_to_oe(r_, v_, mu_)
    %magnitude of r and v
    rmag_ = norm(r_);
    vmag_ = norm(v_);

    %semimajor axis
    a_ = 1/(2/rmag_ - (vmag_^2)/mu_);

    %eccentricity
    evec_ = (vmag_^2/mu_ - 1/rmag_)*r_ - (1/mu_)*dot(r_, v_)*v_;
    e_ = norm(evec_);

    %angular momentum
    hvec_ = cross(r_, v_);
    h_ = norm(hvec_);

    %define orientation vectors
    I_ = [1., 0., 0.];
    J_ = [0., 1., 0.];
    K_ = [0., 0., 1.];

    %node vector
    nvec_ = cross(K_, hvec_/h_);
    n_ = norm(nvec_);

    %inclination
    i_ = acos(dot(hvec_/h_, K_));

    %RAAN
    RAAN_ = acos(dot(nvec_, I_)/n_);
    if dot(nvec_, J_)<0
        RAAN_ = 2*pi - RAAN_;
    end

    %argument of periapse
    w_ = acos(dot(nvec_, evec_)/(n_*e_));
    if dot(evec_, K_)<0
        w_ = 2*pi - w_;
    end

    %true anomaly
    f_ = acos(dot(r_, evec_)/(rmag_*e_));
    if dot(r_, v_)<0
        f_ = 2*pi - f_;
    end

    %convert radians to degrees
    RAAN_ = RAAN_*180/pi;
    w_ = w_*180/pi;
    f_ = f_*180/pi;
    i_ = i_*180/pi;

    a_ = double(a_);
    e_ = double(e_);
    i_ = double(i_);
    RAAN_ = double(RAAN_);
    w_ = double(w_);
    f_ = double(f_);

end