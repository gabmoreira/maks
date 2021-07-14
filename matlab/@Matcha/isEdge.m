function inEdge = isEdge(obj, octave, scale, i, j)

% Compute 2x2 Hessian values from pixel differences 
H00 = obj.dog{octave}(i-1,j,scale) - 2.0 * obj.dog{octave}(i,j,scale) ...
    + obj.dog{octave}(i+1,j,scale);

H11 = obj.dog{octave}(i,j-1,scale) - 2.0 * obj.dog{octave}(i,j,scale) ...
    + obj.dog{octave}(i,j+1,scale);

H01 = ( (obj.dog{octave}(i+1,j+1,scale) - obj.dog{octave}(i+1,j-1,scale)) ...
    - (obj.dog{octave}(i-1,j+1,scale) - obj.dog{octave}(i-1,j-1,scale)) ) / 4.0;

% Compute determinant and trace of the Hessian 
determinant = H00 * H11 - H01 * H01;
trace       = H00 + H11;		


% To detect an edge response, we require the ratio of smallest
% to largest principle curvatures of the DOG function
% (eigenvalues of the Hessian) to be below a threshold.  For
% efficiency, we use Harris' idea of requiring the determinant to
% be above par.EdgeThresh times the squared trace, as for eigenvalues
% A and B, det = AB, trace = A+B.  So if A = 10B, then det = 10B**2,
% and trace**2 = (11B)**2 = 121B**2, so par.EdgeThresh = 10/121 =
% 0.08 to require ratio of eigenvalues less than 10.
	
if (numel(octave) == 1)
    inEdge = ~boolean(determinant > obj.params.edgeThreshold1 * trace * trace);
else
    inEdge = ~boolean(determinant > obj.params.edgeThreshold * trace * trace);
end

end


