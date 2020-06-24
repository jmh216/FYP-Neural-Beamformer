classdef FreeFieldHearingAidArray4 < FreeFieldArray & BinauralArray
        properties (SetAccess=protected)
            % properties from ElobesMicArray
            sensorCartesianPositionsDefault
            refChan
    
            % properties from BinauralArray
            refChanLeft
            refChanRight
            channelsLeft
            channelsRight
        end
    methods
        function[obj] = FreeFieldHearingAidArray4()
       
            % Use superclass to create the object
            obj = obj@FreeFieldArray();
            
            % Populate the parameters
            obj.sensorCartesianPositionsDefault = predefinedSensorPositions();
            obj.refChan = 0;  % reference is the origin
            obj.refChanLeft = 1;
            obj.refChanRight = 2;
            obj.channelsLeft = [1;3];
            obj.channelsRight = [2;4];
        end
    end
    
end

function[sensor_pos] = predefinedSensorPositions()
% evaluate to determine the postitions of the elements relative to
% the origin
radius = 0.1;     % radius on which microphones lie [metres]
spacing = 0.02;   % distance between front and back microphones [metres]

angle_offset = asin(spacing/2 /radius);
az = [pi/2 - angle_offset; ...%front left
    -pi/2 + angle_offset; ...%front right
    pi/2 + angle_offset; ...%rear left
    -pi/2 - angle_offset];   %rear right
inc = pi/2 * ones(4,1);

sensor_pos = radius * [cos(az).*sin(inc), sin(az).*sin(inc), cos(inc)]; % [x,y,z] offsets of sensors
end