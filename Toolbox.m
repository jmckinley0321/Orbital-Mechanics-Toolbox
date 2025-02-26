%% Orbtial Mechanics Toolbox
% Jack McKinley - Purdue Astrodynamics
%%

classdef Toolbox
    methods(Static)

        % Function for converting km to AU
        function AU = km2AU(km)
            % Inputs
                % km         - Distance [km]
            % Outputs
                % AU         - Distance [AU]

            % Conversion Factor
            CF = 149597898;

            % Conversion
            AU = km/CF;
        end

        % Function for converting AU to km
        function km = AU2km(AU)
            % Inputs
                % AU         - Distance [AU]
            % Outputs
                % km         - Distance [km]

            % Conversion Factor
            CF = 149597898;

            % Conversion
            km = AU*CF;
        end

        % Function for converting seconds to Julian Year, Day
        function jy = s2jy(s)
            % Inputs
                % s         - Time [s]
            % Outputs
                % jy        - Time [JY, JD]
    
            % Conversion Factors
            j_dy = 365.25; % Julian Days in a Julian Year
            j_secd = 3600*24; % Seconds in a Julian day
            j_secy = j_secd*j_dy; % Seconds in a Julian Year
        
            % Conversion
            jy = [floor(s/(j_secy)) (mod(s, j_secy)/(j_secd))];
        end

        % Convert r-theta-h to e-p-h frame
        function [eph] = rth2eph(rth, ta)
            % Inputs:
                % rth         - Vector in r-theta-h frame
                % ta          - True anomaly [rad]
            % Outputs
                % eph         - Vector in e-p-h frame
        
            % Convert between r-theta-h and e-p-h frames
            e = rth(1)*cos(ta) - rth(2)*sin(ta) + rth(3)*0;
            p = rth(1)*sin(ta) + rth(2)*cos(ta) + rth(3)*0;
            h = rth(1)*0 + rth(2)*0 + rth(3)*1;
        
            eph = [e p h];
        end
        
        % Convert r-theta-h to x-y-z frame
        function [xyz] = rth2xyz(rth, ta, i, raan, aop)
            % Inputs:
                % rth         - Vector in r-theta-h frame
                % ta          - True anomaly [rad]
                % i           - Inclination [rad]
                % raan        - Right ascension of ascencing node [rad]
                % aop         - Argument of periapsis [rad]
            % Outputs:
                % xyz         - Vector in x-y-z frame 
        
            % Convert true anomaly to theta for frame conversion
            theta = ta + aop;
        
            % Convert between r-theta-h and x-y-z frames
            x = rth(1)*(cos(raan)*cos(theta) - sin(raan)*cos(i)*sin(theta)) ...
                + rth(2)*(-cos(raan)*sin(theta) - sin(raan)*cos(i)*cos(theta)) ...
                + rth(3)*(sin(raan)*sin(i));
            y = rth(1)*(sin(raan)*cos(theta) + cos(raan)*cos(i)*sin(theta)) ...
                + rth(2)*(-sin(raan)*sin(theta) + cos(raan)*cos(i)*cos(theta)) ...
                + rth(3)*(-cos(raan)*sin(i));
            z = rth(1)*(sin(i)*sin(theta)) ...
                + rth(2)*(sin(i)*cos(theta)) ...
                + rth(3)*(cos(i));
        
            xyz = [x y z];
        end

        % Convert eccentric anomaly to true anomaly
        function [ta0] = ea2ta(E0, a, e, p)
            % Inputs:
                % E0          - Eccentric anomaly [rad]
                % a           - Semimajor axis [km]
                % e           - Eccentricity
                % mu          - Semilatus rectum [km]
            % Outputs:
                % ta0         - True anomaly [rad]

            % Calculate radius from eccentric anomaly
            r0 = a*(1 - e*cos(E0));
            disp(r0)

            % Calculate true anomaly from calculated radius
            ta0 = acos(((p/r0) - 1)/e);
            if (E0 > pi) % Quadrant check
                ta0 = -ta0;
            end
        end

        % Solve ellipse given semi major axis and eccentricity
        function [rp, ra, p, h, period, n, epsilon] = solve_ellipse(a, e, mu)
            % Inputs:
                % a           - Semimajor axis [km]
                % e           - Eccentricity
                % mu          - Gravitational parameter [km^3/s^2]
            % Outputs:
                % rp          - Periapsis distance [km]
                % ra          - Apoapsis distance [km]
                % p           - Semilatus rectum [km]
                % h           - Angular momentum [kg*km^2/s]
                % period      - Orbital period [s]
                % n           - Mean motion [rad/s]
                % epsilon     - Orbital energy [km^2/s^2]
            
            % Calculate periapsis distance 
            rp = a*(1 - e);

            % Calculate apoapsis distance
            ra = a*(1 + e);

            % Calculate semilatus rectum
            p = a*(1 - e^2);

            % Calculate angular momentum
            h = sqrt(mu*p);

            % Calculate period
            period = 2*pi*sqrt((a^3)/mu);

            % Calculate mean motion
            n = sqrt(mu/(a^3));

            % Calculate orbital energy
            epsilon = -mu/(2*a);
        end

        % Solve hyperbola given semi major axis and eccentricity
        function [rp, p, b, d, h, ta_inf, v_inf, flyby, epsilon] = solve_hyperbola(a, e, mu)
            % Inputs:
                % a           - Semimajor axis [km]
                % e           - Eccentricity
                % mu          - Gravitational parameter [km^3/s^2]
            % Outputs:
                % rp          - Periapsis distance [km]
                % p           - Semilatus rectum [km]
                % b           - Semi minor axis [km]
                % d           - Aim point [km]
                % h           - Angular momentum [kg*km^2/s]
                % ta_inf      - True anomaly at infinity [rad]
                % v_inf       - Hyperbolic escape velocity [km/s]
                % flyby       - Flyby angle [rad]
                % epsilon     - Orbital energy [km^2/s^2]

            % Calculate periapsis distance
            rp = a*(e - 1);

            % Calculate semilatus rectum
            p = a*((e^2) - 1);

            % Calculate semi minor axis
            b = a*sqrt((e^2) - 1);

            % Calculate aim point
            d = sqrt(((rp + a)^2) - b^2);

            % Calculate angular momentum
            h = sqrt(mu*p);

            % Calculate true anomaly at infinity
            ta_inf = acos(-1/e);

            % Calculate hyperbolic escape velocity
            v_inf = sqrt(mu/a);

            % Calculate flyby angle
            flyby = 2*asin(1/e);

            % Calculate orbital energy
            epsilon = mu/(2*a);
        end

        % Solve state given true anomaly or radius, ascending/descending and conic 
        function [r0, ta0, v0, gamma0, E0, H0, M0, t_since_p] = solve_state(a, e, mu, r0, ta0, asc_desc, conic)
            % Inputs:
                % a           - Semimajor axis [km]
                % e           - Eccentricity
                % mu          - Gravitational parameter [km^3/s^2]
                % r0          - Radius at state [km] (Optional)
                % ta0         - True anomaly at state [rad] (Optional)
                % asc_desc    - Ascending or descending on orbit (Optional)
                % conic       - Type of conic section (Default Ellipse)
            % Outputs:
                % r0          - State orbital radius [km]
                % ta0         - State true anomaly [rad]
                % v0          - State velocity [km/s]
                % gamma0      - State flight path angle [rad]
                % E0          - State eccentric anomaly [rad]
                % H0          - State hyperbolic anomaly
                % M0          - State mean anomaly [rad]
                % ta_inf      - State mean anomaly [rad]
                % t_since_p   - Time since periapsis [s]

            % Calculate semilatus rectum depending on conic section
            if strcmp(conic, "Hyperbola")
                p = a*((e^2) - 1);
            else
                p = a*(1 - (e^2));
            end

            % Calculate angular momentum
            h = sqrt(p*mu);

            % Calculate either ta0 or t0 depending on missing argument
            if isnan(ta0)
                ta0 = acos(((p/r0) - 1)/e);
                if strcmp(asc_desc, "Ascending")
                    ta0 = abs(ta0);
                elseif strcmp(asc_desc, "Descending")
                    ta0 = -abs(ta0);
                end
            elseif isnan(r0)
                r0 = p/(1 + e*cos(ta0));
            end

            % Calculate state velocity
            v0 = sqrt(mu*((2/r0) - (1/a)));

            % Calculate state flight path angle
            gamma0 = acos(h/(r0*v0));
            if strcmp(asc_desc, "Ascending")
                gamma0 = abs(gamma0);
            elseif strcmp(asc_desc, "Descending")
                gamma0 = -abs(gamma0);
            end

            % Calculate state eccentric/hyperbolic anomaly depending on conic           
            if strcmp(conic, "Hyperbola")
                H0 = cosh((a + r0)/(e*a));
                if strcmp(asc_desc, "Descending")
                    H0 = -H0;
                end
                E0 = NaN;
            else
                E0 = acos((cos(ta0) + e)/(1 + e*cos(ta0)));
                if strcmp(asc_desc, "Descending")
                    E0 = (2*pi) - E0;
                end
                H0 = NaN;
            end

            % Calculate mean anomaly if ellipse
            if strcmp(conic, "Hyperbola")
                M0 = NaN;
            else
                M0 = E0 - e*sin(E0);
            end

            % Calculate time since periapsis depending on conic
            if strcmp(conic, "Hyperbola")
                t_since_p = sqrt((a^3)/mu)*(H0 - e*sin(H0));
            else 
                t_since_p = sqrt((a^3)/mu)*(E0 - e*sin(E0));
            end
        end

        % Function to calcualte Hohmann Transfer without Local Gravity Fields
        function [dv_Hohh, dv_dep, dv_arr, TOF, phi, t_s] = hohmann_nograv(periapsis, apoapsis, mu)
            % Inputs:
                % periapsis           - Periapsis of Hohmann Transfer [km]
                % apoapsis            - Apoapsis of Hohmann Transfer [km]
                % mu                  - Gravitational Parameter [km^3/s^2]
            % Outputs:
                % dv_total            - ΔV Total for Hohmann Transfer [km/s]
                % dv_dep              - Departure ΔV [km/s]
                % dv_arr              - Arrival ΔV [km/s]
                % TOF                 - Time of flight [s]
                % phi                 - Phase angle [rad]
                % t_s                 - Synodic period [s]

            % Calculate semi major axis of transfer orbit
            a_Hohh = (periapsis + apoapsis)/2; 

            % Calculate pre-maneuver departure velocity of transfer orbit
            vdep_minus = sqrt(mu/periapsis);

            % Calculate post-maneuver departure velocity of transfer orbit
            vdep_plus = sqrt(mu*((2/periapsis) - (1/a_Hohh)));

            % Calculate departure ΔV
            dv_dep = abs(vdep_plus - vdep_minus);

            % Calculate pre-maneuver arrival velocity of transfer orbit
            varr_minus = sqrt(mu*((2/apoapsis) - (1/a_Hohh)));

            % Calculate post-maneuver arrival velocity of transfer orbit
            varr_plus = sqrt(mu/apoapsis);

            % Calculate arrival ΔV
            dv_arr = abs(varr_plus - varr_minus);

            % Calculate total ΔV
            dv_Hohh = dv_dep + dv_arr;

            % Calculate period of transfer orbit
            period_Hohh = 2*pi*sqrt((a_Hohh^3)/mu);

            % Calculate time of flight
            TOF = period_Hohh/2;

            % Calculate mean motion of pre and post-transfer orbits
            n_minus = sqrt(mu/(periapsis^3));
            n_plus = sqrt(mu/(apoapsis^3));

            % Calculate phase angle
            phi = pi - (n_plus*TOF);

            % Calculate synodic period
            t_s = (2*pi)/(n_minus - n_plus);
        end

        function [a, alpha, beta] = lambert_solver(s, c, TOF, mu, type)
            % Inputs:
                % s                   - Semiperimeter of Space Triangle [km]
                % c                   - Chord of Space Triangle [km]
                % TOF                 - Specified time of flight [s]
                % mu                  - Gravitational Parameter [km^3/s^2]
                % type                - Transfer Type
            % Outputs:
                % a                   - Transfer semimajor axis [km]
                % alpha               - Defined quantity for Lambert [rad]
                % beta                - Defined quantity for Lambert [rad]

            % Calculate alpha0, beta0 for transfer type
            if (strcmp(type, "1H") || strcmp(type, "2H"))
                falpha0 = @(a) 2*asinh(sqrt((s/(2*a))));
                fbeta0 = @(a) 2*asinh(sqrt(((s - c)/(2*a))));
            else
                falpha0 = @(a) 2*asin(sqrt((s/(2*a))));
                fbeta0 = @(a) 2*asin(sqrt(((s - c)/(2*a))));
            end
    
            % Correct alpha, beta for transfer type
            if (strcmp(type, "1A") || strcmp(type, "1H"))
                falpha = @(a) falpha0(a);
                fbeta = @(a) fbeta0(a);
            elseif (strcmp(type, "1B"))
                falpha = @(a) 2*pi - falpha0(a);
                fbeta = @(a)  fbeta0(a);
            elseif (strcmp(type, "2A") || strcmp(type, "2H"))
                falpha = @(a) falpha0(a);
                fbeta = @(a) -fbeta0(a);
            elseif (strcmp(type, "2B"))
                falpha = @(a) (2*pi) - falpha0(a);
                fbeta = @(a) -fbeta0(a);
            end 
    
            % Build hyperbolic and elliptical lambert equations
            lambert_hyp = @(a) (sqrt(mu)*TOF) - a^(3/2)*((sinh(falpha(a)) - falpha(a)) - (sinh(fbeta(a)) - fbeta(a)));
            lambert_ell = @(a) (sqrt(mu)*TOF) - a^(3/2)*((falpha(a) - fbeta(a)) - (sin(falpha(a)) - sin(fbeta(a))));

            % Initialize guess for semimajor axis
            a_guess = 2*s;

            % Numerically solve for semimajor axis by transfer type
            if (strcmp(type, "1H") || strcmp(type, "2H"))
                a = real([fsolve(lambert_hyp, a_guess)]);
            else
                a = real([fsolve(lambert_ell, a_guess)]);
            end

            % Set values for alpha, beta 
            alpha = falpha(a);
            beta = fbeta(a);
        end
    end
end
