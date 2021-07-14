classdef Posegraph < handle

    properties (Access = public)
        edges;
        nodes;
        
        numEdges;
        numNodes;
    end
    
    properties (Access = private)
        filename;
        nodeFields;
        edgeFields;
    end
    
    methods
        % Class constructor
        function obj = Posegraph()
            obj.numEdges   = 0;
            obj.numNodes   = 0;
            obj.nodeFields = "id";
            obj.edgeFields = "id";
        end
        
        % Reset graph
        clear(obj);
        
        % Specify node fields
        configNodes(obj, varargin);
   
        % Check stuff
        idx = containsEdge(obj, id);
        idx = containsNode(obj, id);
        
        % Graph changing stuff
        addEdge(obj, id, varargin);
        addNode(obj, id, varargin);

        % Nodewise operation
        nodeOp(obj, op, field);
                
        % Algebraic graph structure
        A = adjacency(obj, varargin);    
        D = degree(obj, varargin);       
        L = laplacian(obj, varargin);
        f = fiedler(obj);   
    end

end

