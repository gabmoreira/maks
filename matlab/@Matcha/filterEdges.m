function filterEdges(obj)
% FILTEREDGES Eliminate keypoints with a high edge response
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

fprintf("Matcha: Computing edge response - [")
countRemovals = 0;
removeIdx = [];

consoleCounter = 0;
for k=1:size(obj.extrema, 1)
    if (round(100*k/size(obj.extrema, 1)) > consoleCounter)
        fprintf("=");
        consoleCounter = consoleCounter + 10;
    end
    
    if( obj.isEdge(obj.extrema(k,1), obj.extrema(k,2), obj.extrema(k,3), obj.extrema(k,4)) )
        removeIdx = [removeIdx; k];
        countRemovals = countRemovals + 1;
    end
end

if (~isempty(removeIdx))
    obj.extrema(removeIdx,:) = [];
end
        
fprintf("] - Eliminated %d edges. \n", countRemovals);
end

