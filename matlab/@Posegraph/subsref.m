function varargout = subsref(obj, s)
% SUBSREF Overload of MATLAB's internal subscript functions
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

switch s(1).type
    % Let MATLAB handle this one
    case '.'
        [varargout{1:nargout}] = builtin('subsref',obj,s);
        
    % Node / Edge query
    case '()'
        % Handle nodes
        if length(s(1).subs) == 1
            
            % Query node id
            id  = s(1).subs{1};
            % Look for the index of this id in the graph
            idx = obj.containsNode(id);
            
            if (idx <= obj.numNodes)
                node = {};
                for field=obj.nodeFields
                    node.(field) = obj.nodes.(field){idx};
                end
                varargout{1} = node;
                
                if length(s) > 1
                    [varargout{1:nargout}] = builtin('subsref',varargout{1},s(2:end));
                end
            else
                error('Node not in the graph');
            end

        % Handle edges
        elseif length(s(1).subs) == 2 
            i = s(1).subs{1};
            j = s(1).subs{2};

            idx = obj.containsEdge([i,j]);
            if (idx <= obj.numEdges)
                edge = {};
                for field=obj.nodeFields
                    edge.(field) = obj.edges.(field){idx};
                end
                varargout{1} = edge;
                
                if length(s) > 1
                    [varargout{1:nargout}] = builtin('subsref',varargout{1},s(2:end));
                end
            else
                error('Edge not in the graph');
            end
        else
            [varargout{1:nargout}] = builtin('subsref',obj,s);
        end
        
    % Let MATLAB handle this one
    case '{}'
        [varargout{1:nargout}] = builtin('subsref',obj,s);
        
    % Don't even try
    otherwise
     error('Not a valid indexing expression')
end

end