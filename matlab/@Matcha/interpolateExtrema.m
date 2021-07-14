function interpolateExtrema(obj)
% INTERPOLATEEXTREMA Scale-space interpolation of all DoG extrema
%
% Other m-files required: Matcha.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

assert(obj.COMPUTED_DOG_EXTREMA, 'Matcha: Must compute Difference-of-Gaussian (dog) extrema first.');

fprintf("Matcha: Interpolating extrema   - [");

% Some extrema will be removed if the interpolation goes bad
removeIdx = [];

% Allocate space for x, y, sigma, peakval
obj.extrema = [obj.extrema, zeros(size(obj.extrema,1),4)];

consoleCounter = 0;
for k=1:size(obj.extrema, 1)
    if (round(100*k/size(obj.extrema, 1)) > consoleCounter)
        fprintf("=");
        consoleCounter = consoleCounter + 10;
    end
    
    
    [success, s, i, j, x, y, sigma, peakval] = obj.interpolatePoint(obj.extrema(k,1), ...
                                                                    obj.extrema(k,2), ...
                                                                    obj.extrema(k,3), ...
                                                                    obj.extrema(k,4));
    if (success)
        obj.extrema(k,2:8) = [s, i, j, x, y, sigma, peakval];
    else
        removeIdx = [removeIdx; k];
    end

end

obj.extrema(removeIdx,:) = [];
fprintf("] - %d interpolated extrema.\n", size(obj.extrema, 1));
end

