function writeG2O(filename, M, A, varargin)
% WRITEG2O Save data as .g2o file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

valid_filename = @(x) isstring(x);
valid_M        = @(x) issparse(x);
valid_A        = @(x) issparse(x);
valid_group    = @(x) isstring(x);
valid_cycle    = @(x) islogical(x);

default_cycle = false;
default_group = "SE3";

p = inputParser;

addRequired(p, 'filename', valid_filename);
addRequired(p, 'M', valid_M);
addRequired(p, 'A', valid_A);

addParameter(p, 'group', default_group, valid_group);
addParameter(p, 'cycle', default_cycle, valid_cycle);

parse(p, filename, M, A, varargin{:});

fileID = fopen(filename + ".g2o",'w');

n = size(A,1);

% If graph has one cycle, saves the edges as 1-2, 2-3, 3-4, ..., n-1
if (p.Results.cycle)
    fprintf("Cycle flag on.\n")
    for k = 1:n-1       
        switch p.Results.group
            case 'SE3'
                pose = full(M(k*4-3:k*4,(k+1)*4-3:(k+1)*4));
                rot  = pose(1:3,1:3);
                t    = pose(1:3,4);
            case 'SO3'
                rot = full(M(k*3-2:k*3,(k+1)*3-2:(k+1)*3));
                t   = zeros(3,1);
        end
        
        q  = rotToQuat(rot);
        qx = round(q(1),15);
        qy = round(q(2),15);
        qz = round(q(3),15);
        qw = round(q(4),15);
        
        fprintf(fileID,'EDGE_SE3:QUAT %d %d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', ...
            k-1, k, t(1), t(2), t(3), qx, qy, qz, qw, ...
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
    end

    switch p.Results.group
    	case 'SE3'
            pose = full(M(end-3:end,1:4));
            rot  = pose(1:3,1:3);
        	t    = pose(1:3,4);
        case 'SO3'
            rot = full(M(end-2:end,1:3));
        	t   = zeros(3,1);
    end
    q  = rotToQuat(rot);
    qx = round(q(1),15);
    qy = round(q(2),15);
    qz = round(q(3),15);
    qw = round(q(4),15);
    fprintf(fileID,'EDGE_SE3:QUAT %d %d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', ...
            n-1, 0, t(1), t(2), t(3), qx, qy, qz, qw, ...
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);

% Graph is not a cycle       
else
    A = triu(A,1);
    [ii,jj] = find(A);
     fprintf("Writing %d edges\n", numel(ii));

    for k = 1:numel(ii)
        i = ii(k);
        j = jj(k);
        
        switch p.Results.group
            case 'SE3'
                pose = full(M(i*4-3:i*4,j*4-3:j*4));
                rot  = pose(1:3,1:3);
                t    = pose(1:3,4);
            case 'SO3'
                rot = full(M(i*3-2:i*3,j*3-2:j*3));
                t   = zeros(3,1);
        end

        q  = rotToQuat(rot);
        qx = round(q(1),15);
        qy = round(q(2),15);
        qz = round(q(3),15);
        qw = round(q(4),15);

        fprintf(fileID,'EDGE_SE3:QUAT %d %d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', ...
            i-1, j-1, t(1), t(2), t(3), qx, qy, qz, qw, ...
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
    end
end

fclose(fileID);
fprintf("File saved as %s\n", filename + ".g2o");

end
