function [offset, peakval, info] = fitQuadratic(obj, o, s, i, j)
% Apply the method developed by Matthew Brown (see BMVC 02 paper) to
% fit a 3D quadratic function through the DOG function values around
% the location (o,s,i,j) at which a peak has
% been detected.  Return the interpolated peak position as a vector
% in "offset", which gives offset from position (s,i,j).  The
% returned value is the interpolated DOG magnitude at this peak.
info = true;

% Fill in the values of the gradient from pixel differences
g = zeros(3,1);
g(1) = (obj.dog{o}(i,j,s+1) - obj.dog{o}(i,j,s-1)) / 2.0;
g(2) = (obj.dog{o}(i,j+1,s) - obj.dog{o}(i,j-1,s)) / 2.0;
g(3) = (obj.dog{o}(i+1,j,s) - obj.dog{o}(i-1,j,s)) / 2.0;

% Fill in the values of the Hessian from pixel differences 
H = zeros(3,3);

H(1,1) = obj.dog{o}(i,j,s-1) + obj.dog{o}(i,j,s+1) - 2.0 * obj.dog{o}(i,j,s);
H(2,2) = obj.dog{o}(i,j-1,s) + obj.dog{o}(i,j+1,s) - 2.0 * obj.dog{o}(i,j,s);
H(3,3) = obj.dog{o}(i-1,j,s) + obj.dog{o}(i+1,j,s) - 2.0 * obj.dog{o}(i,j,s);

H(1,2) = ( obj.dog{o}(i,j+1,s+1) - obj.dog{o}(i,j-1,s+1) - obj.dog{o}(i,j+1,s-1) + obj.dog{o}(i,j-1,s-1) ) / 4.0;
H(1,3) = ( obj.dog{o}(i+1,j,s+1) - obj.dog{o}(i-1,j,s+1) - obj.dog{o}(i+1,j,s-1) + obj.dog{o}(i-1,j,s-1) ) / 4.0;
H(2,3) = ( obj.dog{o}(i+1,j+1,s) - obj.dog{o}(i-1,j+1,s) - obj.dog{o}(i+1,j-1,s) + obj.dog{o}(i-1,j-1,s) ) / 4.0;

H(2,1) = H(1,2);
H(3,2) = H(2,3);
H(3,1) = H(1,3);

% Solve the 3x3 linear sytem, Hx = -g. Result, x, gives peak offset.
offset = - (H \ g);

if ( any(isnan(offset(:))) )
    info = false;
end

% return value of DOG at peak location using initial value plus
peakval = obj.dog{o}(i,j,s) - 0.5 * offset' * g;
end
