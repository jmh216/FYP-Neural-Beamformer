% https://uk.mathworks.com/help/phased/ref/mvdrweights.html
% Phased Array System Toolbox

N = 6;
d = 0.05;
elementPos = (0:N-1)*d;
% Sn  = sensorcov(elementPos,[0 0],db2pow(-10));
Sn  = sensorcov(elementPos,[30 45],db2pow(-10));

% w = mvdrweights(elementPos,[30 45],Sn);
w = mvdrweights(elementPos,0,Sn);

plotangl = -90:90;
vv = steervec(elementPos,plotangl);
plot(plotangl,mag2db(abs(w'*vv)))
xlim([-90 90]);
grid on
xlabel('Azimuth Angle (degrees)');
ylabel('Normalized Power (dB)');
legend('30 deg','45 deg');
title('MVDR Array Pattern')