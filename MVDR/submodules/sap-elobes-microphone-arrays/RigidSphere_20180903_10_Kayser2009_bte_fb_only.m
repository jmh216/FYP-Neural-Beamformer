classdef RigidSphere_20180903_10_Kayser2009_bte_fb_only < RigidSphereArray & BinauralArray
    properties (SetAccess=protected)
        % properties from ElobesMicArray already implemented in SphericalHarmonicSoundFieldArray
        
        % properties from BinauralArray
        refChanLeft
        refChanRight
        channelsLeft
        channelsRight
        
        % new properties to allow level to be calibrated to measurements
        micGain
    end
    methods
        function[obj] = RigidSphere_20180903_10_Kayser2009_bte_fb_only()
            
            % use superclass to create the object with input parameter
            obj = obj@RigidSphereArray(0.09);
            
            % override the default properties of the superclass
            obj.sensorCartesianPositionsDefault = predefinedSensorPositions();
            obj.refChan = 0; % reference is the origin
            
            % populate the parameters
            obj.refChanLeft = 1;
            obj.refChanRight = 2;
            obj.channelsLeft = [1;3];
            obj.channelsRight = [2;4];
            
            obj.micGain = 10^(-32/20) * ones(1,size(obj.sensorCartesianPositionsDefault,1)); %[1, nChans]
        end
        
        function[H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin)
            % use supercalss to do the work
            [H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc@RigidSphereArray(obj,src_az,src_inc,varargin);
            
            % apply microphone gain
            H = bsxfun(@times,obj.micGain,H);
                        
        end
    end

end

function[sensor_pos] = predefinedSensorPositions()
% evaluate to determine the postitions of the elements relative to
% the origin
radius = 0.09;%obj.sphereRadius; %0.082;     % radius on which microphones lie [metres]
spacing = 0.015;   % distance between front and back microphones [metres]
angle_offset = deg2rad(7.5); % line bisecting microphones is offset back from a line through the sphere

angle_spacing = asin(spacing/2 /radius);
az = [pi/2 + angle_offset - angle_spacing; ...%front left
    -(pi/2 + angle_offset) + angle_spacing; ...%front right
    pi/2 + angle_offset + angle_spacing; ...%rear left
    -(pi/2 + angle_offset) - angle_spacing];   %rear right
inc = pi/2 * ones(4,1);

sensor_pos = radius * [cos(az).*sin(inc), sin(az).*sin(inc), cos(inc)]; % [x,y,z] offsets of sensors
end