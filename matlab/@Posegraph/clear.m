function clear(obj)
% CLEAR Resets graph
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

obj.numEdges = 0;
obj.numNodes = 0;

obj.nodeFields = "id";
obj.edgeFields = "id";

end

