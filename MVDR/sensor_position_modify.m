% classdef FreeField_20180202_01_ULA_configurable < FreeFieldArray
%         properties (SetAccess=protected)
%             % properties from ElobesMicArray
%             sensorCartesianPositionsDefault
%             refChan
%     
%             % properties from BinauralArray
%             refChanLeft
%             refChanRight
%             channelsLeft
%             channelsRight
%         end
%     methods
%         function[obj] = FreeField_20180202_01_ULA_configurable(nSensors,sensorSpacing)
%        
%             % Use superclass to create the object
%             obj = obj@FreeFieldArray();
%             
%             % Populate the parameters
%             obj.sensorCartesianPositionsDefault = predefinedSensorPositions(nSensors,sensorSpacing);
%             obj.refChan = 1;  % reference is the origin
%         end
%     end
%     
% end

nSensors = 6;
sensorSpacing = 0.10;

for i = 1:3
x_pos = 0.10 * (0:-1:-(3-1)).';
y_pos = zeros(3,1);
z_pos = zeros(3,1);
end
sensor_pos = [x_pos, y_pos, z_pos ];
for i = 4:6
x_pos_bot = sensorSpacing * (0:-1:-(3-1)).';
y_pos_bot = zeros(3,1);
z_pos_bot = 0.19 * ones(3,1);
end
sensor_pos_bot = [x_pos_bot,y_pos_bot,z_pos_bot]; % [x,y,z] offsets of sensors
final_pos = [sensor_pos; sensor_pos_bot];

figure;
plot3(final_pos(:,1),final_pos(:,2),final_pos(:,3),'o');
hold on;
Room = [8 8 4];			% Room dimensions [x y z] (m) 
axis([-1 1 -1 1 -1 1]);
xlabel("x-axis/m"); ylabel("y-axis/m"); zlabel("z-axis/m");
legend('Microphone')
title("Position of microphone array");
grid on; box on; axis square;hold off;

% function[sensor_pos] = predefinedSensorPositions(nSensors,sensorSpacing)
% % evaluate to determine the postitions of the elements relative to
% % the origin
% % to number from left to right from listener perspective align with y-axis
% 
% %halfN = (nSensors-1)/2;
% %if rem(halfN,1)~=0, error('nSensors must be odd'),end
% x_pos = sensorSpacing * (0:-1:-(nSensors-1)).';
% y_pos = zeros(nSensors,1);
% z_pos = zeros(nSensors,1);
% 
% sensor_pos = [x_pos,y_pos,z_pos]; % [x,y,z] offsets of sensors
% end