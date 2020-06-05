classdef WaveGuide < handle
    %WAVEGUIDE Class for encapuslating a laser waveguide
    % 
    % This class encapsulates a laser waveguide in terms of the core and
    % cladding refractive indexes, width (or radius) of the guide and the
    % free-space wavevector of the laser.
    %
    % On initialisation, the propagation constant (beta) and effective
    % index (neff) are calculated for the lowest supported mode of the
    % guide.
    %
    % Presently, this functionality is only available for slab guides.
    %
      
    properties
        
        ncore	{mustBeNumeric} % Core refractive index
        nclad	{mustBeNumeric} % Cladding refractive index
        width	{mustBeNumeric} % Width of slab waveguide (micron)
        radius	{mustBeNumeric} % Radius of cylindrical waveguide (micron)
        k0      {mustBeNumeric} % Free-space wavevector (1/micron)
        Gamma   {mustBeNumeric} % Optical confinement factor
        type    {mustBeNumeric} % Type (0 = slab; 1 = cylinder)
        
    end
    
    properties (SetAccess=protected)
        
        beta    {mustBeNumeric} % Propagation constant (1/micron)
        neff    {mustBeNumeric} % Effective refractive index
        is_calc = false    % Have dependent properties been calculated? 
        
    end
    
    methods
        
        function obj = WaveGuide(n1, n2, w, k0, varargin)
            %WAVEGUIDE Constructs an instance of the class WaveGuide
            %
            % Usage:
            %
            %	obj = WaveGuide(n1, n2, w, k0[, type]);   
            %
            % Arguments:
            %
            %   n1      Core refractive index
            %   n2      Cladding refractive index
            %   w       Width / Radius (depending on type)
            %   k0      Free-space wavevector (1/micron)
            %   type    [Optional] Sets the type of the waveguide:
            %
            %       0   Slab waveguide (currently only type supported)
            %
            % The constructor also calculates and sets the dependent
            % properties:
            %
            %   beta    Propagation constant (1/micron)
            %   neff    Effective refractive index
            %   Gamma   Optical confinement factor
            %
            % An internal property 'is_calc' is then set to indicate that
            % this calculation has been done.
            
            obj.ncore = n1;
            obj.nclad = n2;
            obj.k0 = k0;
            
            if (nargin == 4)
                % type not passed - set as slab (default)
                obj.width = w;
                obj.radius = 0.0;
                obj.type = 0;
                obj.setSlab();
                
            else
                
                type = varargin{1};
                
                if (type == 0)
                    % Waveguide is a slab
                    obj.width = w;
                    obj.radius = 0.0;
                    obj.type = 0;
                    obj.setSlab();
                else
                    % Not currently supported
                    errID = 'WaveGuide:Exception';
                    errMsg = 'Only slab waveguides are currently supported';
                    ME = MException(errID, errMsg);
                    throw(ME);
                    
                end
            end
        end
        
        function bool = isCalc(obj)
            %ISCALC Tests whether dependent properties have been calculated
            %with current independent properties
            
            bool = obj.is_calc;
            
        end
        
        function reCalc(obj)
            %RECALC Recalculates dependent properties (if necessary)
            
            if (~obj.is_calc)
                
                obj.setSlab();
                
            end
            
        end
        
        %% Setters
        
        function set.ncore(obj, n1)
            %SET.NCORE Sets the core refractive index
            %
            % Usage:
            %
            %   obj.ncore = n1;
            %
            % The setter also sets 'is_calc' to false is ncore has been
            % changed.
            
            if (n1 ~= obj.ncore)
                obj.ncore = n1;
                obj.is_calc = false; %#ok<MCSUP>
            end
            
        end
        
        function set.nclad(obj, n2)
            %SET.NCLAD Sets the cladding refractive index
            %
            % Usage:
            %
            %   obj.nclad = n2;
            %
            % The setter also sets 'is_calc' to false is ncore has been
            % changed.
            
            if (n2 ~= obj.nclad)
                obj.nclad = n2;
                obj.is_calc = false; %#ok<MCSUP>
            end
            
        end
        
        function set.width(obj, w)
            %SET.WIDTH Sets the waveguide width (microns)
            %
            % Usage:
            %
            %   obj.width = w;
            %
            % The setter also sets 'is_calc' to false is ncore has been
            % changed.
            
            if (w ~= obj.width)
                obj.width = w;
                obj.is_calc = false; %#ok<MCSUP>
            end
            
        end
        
        function set.radius(obj, r) %#ok<INUSD>
            %SET.RADIUS Sets the waveguide radius (microns)
            %
            % NOT CURRENTLY SUPPORTED.
            %
            errID = 'WaveGuide:Exception';
            errMsg = 'Only slab waveguides are currently supported';
            ME = MException(errID, errMsg);
            throw(ME);
            
        end
        
        function set.k0(obj, k0)
            %SET.K0 Sets the free-space wave-vector (1/microns)
            %
            % Usage:
            %
            %   obj.k0 = k0;
            %
            % The setter also sets 'is_calc' to false is ncore has been
            % changed.
            
            if (k0 ~= obj.k0)
                obj.k0 = k0;
                obj.is_calc = false; %#ok<MCSUP>
            end
            
        end
        
        function set.Gamma(obj, G)
            %SET.GAMMA Sets the optical confinement factor
            %
            % Usage:
            %
            %   obj.Gamma = G;
            %
            % The setter also sets 'is_calc' to false is ncore has been
            % changed.
            
            if (G ~= obj.Gamma)
                obj.Gamma = G;
                obj.is_calc = false; %#ok<MCSUP>
            end
            
        end
        
        function set.type(obj, type)
            %SET.TYPE Sets the waveguide type
            %
            % Currently redundant
            
            if (type ~= 0)
                
                % Not currently supported
                errID = 'WaveGuide:Exception';
                errMsg = 'Only slab waveguides are currently supported';
                ME = MException(errID, errMsg);
                throw(ME);
                
            else
            
                obj.type = type;
                
            end
             
        end
        
    end
    
    %% Protected methods
    
    methods (Access=protected)
        
        function setSlab(obj) 
            %SETSLAB Calculates the dependent properties for a slab
            %waveguide.
            %
            % These properties are:
            %
            %   beta    Propagation constant (1/micron)
            %   neff    Effective refractive index
            %   Gamma   Optical confinement factor
            %
            % An internal property 'is_calc' is then set to indicate that
            % this calculation has been done.
            
            [b, G] = slabGuideBeta(obj.k0, obj.ncore, obj.nclad, obj.width);
            obj.beta = b;
            obj.Gamma = G;
            obj.neff = obj.beta/obj.k0;
            obj.is_calc = true;
            
        end
    end
end

