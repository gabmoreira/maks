function R = rotrnd(n, sigma)
% ROTRND Random rotations
%
% Axis sampled uniformly over the unit sphere and and angle sampled from
% Normal distribution N(0,sigma).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         g.antunes.moreira (at) gmail.com
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

R = zeros(n*3,3);

rvals     = 2*rand(n,1)-1;
elevation = asin(rvals);
azimuth   = 2*pi*rand(n,1);
[x,y,z]   = sph2cart(azimuth,elevation,1);
magnitude = sigma * randn(n,1);

for i=1:n
    R(i*3-2:i*3,:) = angleAxisToRot([x(i), y(i), z(i)], magnitude(i));
end

if (n > 1)
    R = R * R(1:3,1:3)';
end

end

