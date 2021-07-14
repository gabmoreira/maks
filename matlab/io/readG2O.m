function [e, R, t] = readG2O(filename)
% READG2O Reads edges from G2O file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

fid = fopen(filename, 'r');

R = zeros(3, 3*30000);
t = zeros(3, 30000);
e = zeros(2, 30000);

k = 1;
tline = fgetl(fid);

while ischar(tline)
    tokens = split(tline);
    
    if (strcmp(tokens{1},'EDGE_SE3:QUAT'))
        
        new_edge_ij = [str2double(tokens{2}); str2double(tokens{3})];
        new_edge_ji = [str2double(tokens{3}); str2double(tokens{2})];

        if ( (~ismember(new_edge_ji.', e.', 'rows')) && ...
             (~ismember(new_edge_ij.', e.', 'rows')) )
         
            % Store edge
            e(:,k) = new_edge_ij;
            
            % Store translation
            t(:,k) = [str2double(tokens{4});
                      str2double(tokens{5});
                      str2double(tokens{6})];
            
            % Store rotation
            R(:,k*3-2:k*3) = quatToRot([str2double(tokens{7});
                                        str2double(tokens{8});
                                        str2double(tokens{9});
                                        str2double(tokens{10})]);
            k = k + 1;
        end
    end
    
    tline = fgetl(fid);
    
end

fclose(fid);

% Trim and offset indices by 1 since g2o is indexed at 0
e = e(:,1:k-1) + 1;
R = R(:,1:3*k-3);
t = t(:,1:k-1);

nodes = union(e(1,:), e(2,:)); 

numNodes = numel(nodes);
numEdges = size(e,2);

% Convert to cell array for compatibility with pose graph
R = mat2cell(R, 3, ones(1, numEdges)*3);
t = mat2cell(t, 3, ones(1, numEdges));
e = mat2cell(e, 2, ones(1, numEdges));

fprintf("Loaded %d edges corresponding to %d nodes.\n", numEdges, numNodes);

end

