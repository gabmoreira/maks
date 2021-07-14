function angular = chordalToAngular(chordal)
% CHORDALTOANGULAR Rotation chordal distance to angular distance (deg)
%
% Syntax:  [angular] = chordal2angular(chordal)
%
% Inputs:
%    chordal - (Double) Chordal distance ||R1-R2||_F
%
% Outputs:
%    angular - (Double) Angular distance in degrees between R1 and R2
%
% Example: 
%          
%    chordal = 2
%
%    angular = 90
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       http://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

angular = rad2deg(2 * real(asin(chordal / (sqrt(2) * 2))) );

end
