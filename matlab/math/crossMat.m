function dst = crossMat(x)
% CROSSMAT Cross product matrix representation of a vector x
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

dst = [    0, -x(3),  x(2);
        x(3),     0, -x(1);
       -x(2),  x(1),    0];
   
end

