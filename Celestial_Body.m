%% Solar System Definitions
% Jack McKinley - Purdue AAE Fall 2024
%%

classdef Celestial_Body
    properties
        Name                 % Celestial Body Name        
        Radius = []          % Mean Equatorial Radius [km]  
        Mu = []              % Gravitational Parameter [km^3/s^2]
        SMA = []             % Semimajor Axis of Orbit [km]
        ECC = []             % Eccentricity of Orbit 
        INC = []             % Inclination of Orbit to Elliptic [deg]
        Period = []          % Orbital Period [s]
    end

    methods
        % Constructor
        function obj = Celestial_Body(name, radius, mu, sma, ecc, inc, period)
            obj.Name = name;
            obj.Radius = radius;
            obj.Mu = mu;
            obj.SMA = sma;
            obj.ECC = ecc;
            obj.INC = inc;
            obj.Period = period;
        end
    end

    methods (Static)
        % Initialize Celestial Bodies
        function celestial_bodies = predefined()
            celestial_bodies.Sun = Celestial_Body('Sun', 695990, ...
                132712440017.99, [], [], [], []);
            celestial_bodies.Mercury = Celestial_Body('Mercury', 2439.7, ...
                22032.080486418, 57909101, 7600537, 0.20563661, 7.004097902);
            celestial_bodies.Venus = Celestial_Body('Venus', 6051.9, ...
                324858.59882646, 108207284, 19413722, 0.00676399, 3.39465605);
            celestial_bodies.Earth = Celestial_Body('Earth',  6378.1363, ...
                398600.4415, 149597898, 31558205, 0.01673163, 0.00001531);
            celestial_bodies.Luna = Celestial_Body('Luna', 1738.2, ...
                4902.8005821478, 384400, 2360592, 0.0554, 5.16);
            celestial_bodies.Mars = Celestial_Body('Mars', 3397, ...
                42828.314258067, 227944135, 59356281, 0.09336511, 1.84969142);
            celestial_bodies.Jupiter = Celestial_Body('Jupiter', 71492, ...
                126712767.8578, 778279959, 374479305, 0.04853590, 1.30439695);
            celestial_bodies.Saturn = Celestial_Body('Saturn',  60268, ...
                37940626.061137, 1427387908, 930115906, 0.05550825, 2.48599187);
            celestial_bodies.Uranus = Celestial_Body('Uranus', 25559, ...
                5794549.0070719, 2870480873, 2652503938, 0.04685740, 0.77263783);
            celestial_bodies.Neptune = Celestial_Body('Neptune', 25269, ...
                6836534.0638793, 4498337290, 5203578080, 0.00895439, 1.77004347);
            celestial_bodies.Pluto = Celestial_Body('Pluto', 1162, ...
                981.600887707, 5907150229, 7830528509, 0.24885238, 17.14001206);
        end
        
        % Get Celestial Body
        function celestial_body = get(name)
            celestial_bodies = Celestial_Body.predefined();
            
            celestial_body = celestial_bodies.(name);
        end
    end
end