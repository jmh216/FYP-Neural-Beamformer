classdef Sampled_20181119_Kayser2009_bte_fb_only < SampledArray & BinauralArray
    properties (SetAccess=protected)
        % properties from ElobesMicArray
        sensorCartesianPositionsDefault
        refChan
            
        % properties from BinauralArray
        refChanLeft
        refChanRight
        channelsLeft
        channelsRight
        
        % new properties to allow level to be calibrated to measurements
        %micGain
    end
    
    properties (Hidden)
        % Store some precomputed values
        distance = [];
    end
    methods
        function[obj] = Sampled_20181119_Kayser2009_bte_fb_only(distance)
            
            % use superclass to create the object with input parameter
            obj = obj@SampledArray();
            
            obj.supportsRotation=1; % using anechoic measurements 
            obj.availableInterpolationMethods = {'none','nearest_neighbour'};
            obj.interpolationMethod = 'nearest_neighbour';
            
            obj.sensorCartesianPositionsDefault = predefinedSensorPositions();
            obj.refChan = nan; % reference is the origin
            
            % populate the parameters
            obj.refChanLeft = 1;
            obj.refChanRight = 2;
            obj.channelsLeft = [1;3];
            obj.channelsRight = [2;4];
            
            % specific to this array, have a choice of measurment distances
            if nargin==0
                distance = 300;
            end
            obj.distance = distance/100; % in metres
            
        end
        function[rMin, rMax] = getValidSrcRadiusRange(obj)
            rMin = obj.distance;
            rMax = obj.distance;
        end
        function[ir,src_pos,fs] = loadSampledData(obj)
            % TODO: Establish mechanism for defining paths to binary data
%            wav_path = '~/local_databases/Impulse responses/Oldenburg/HRIR_database_wav/hrir/anechoic';
            wav_path = '/Users/amoore1/data/_third_party/Kayser2009/HRIR_database_wav/hrir/anechoic';
            
            
            % define grid in database co-ordinates
            az_deg_vec = [-180:5:175];
            el_deg_vec = [-10:10:20];
            [el_deg_grid,az_deg_grid] = ndgrid(el_deg_vec,az_deg_vec);
            el_deg = el_deg_grid(:);
            az_deg = az_deg_grid(:);
            
            nPos = size(az_deg,1);
            for ipos = nPos:-1:1
                [ir(:,:,ipos),fs] = readwav(fullfile(wav_path,...
                    sprintf('anechoic_distcm_%d_el_%d_az_%d.wav',...
                    obj.distance*100,el_deg(ipos),az_deg(ipos))));
            end
            % remap and select the required channels
            ir = ir(:,[3 4 7 8],:);
            % define directions in cartesian coordinates according to
            % elobes convention
            src_pos = mysph2cart(deg2rad(-az_deg),deg2rad(90-el_deg),obj.distance*ones(nPos,1));     
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