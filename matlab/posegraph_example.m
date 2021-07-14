% POSEGRAPH_EXAMPLE Mini tutorial on the Posegraph class
%
% Other m-files required: ./@Posegraph/*, ./geo/*, ./math/*, ./sparse/*
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

clc;
close all;
clear all;

addpath('geo');
addpath('math');

g = Posegraph();

% Specify the attributes of each node (type/shape agnostic)
g.configNodes("R", "t");

% Specify the attributes of each edge (type/shape agnostic)
g.configEdges("R", "t");

% Add nodes by specifying an id and the attributes specified earlier
g.addNode(1, 'R', rand(3,3), 't', rand(2,1))
g.addNode(2, 'R', rand(3,3), 't', rand(2,1))
g.addNode(3, 'R', rand(3,3), 't', rand(2,1))

% Add edges with attributes specified earlier between existing node ids
g.addEdge([1;2], 'R', rotrnd(1,0.5), 't', rand(3,1));
g.addEdge([1;3], 'R', rotrnd(1,0.5), 't', rand(3,1));

%% Accessing nodes and edges structs

% Node with id 1
n1 = g(1);

% Edge from id 1 to id 2
e12 = g(1,2);

% Error!
g(4,10);

%% Accessing nodes and edges attributes

% Attribute R of node 1. Same thing as n1.R
R1 = g(1).R;

% Attribute R of edge 12. Same thing as e12.R
R12 = g(1,2).R;

% Error ! Non-existing attribute
g(1).y

%% Apply a function to a certain attribute in ALL the nodes 

% Go over all nodes and apply function ./math/projectToSO3 to attribute R
g.nodeOp(@projectToSO3, "R");

%% Graph structure

% Graph degree as a vector
degreeVec = g.degree();

% Graph degree as a sparse diagonal matrix
degreeMat = g.degree('asmatrix', true);

% Sparse graph adjacency matrix
A = g.adjacency();

% Sparse graph adjacency with weights from the "R" atribute of the nodes
AR = g.adjacency('field' , "R");

% Sparse graph Laplacian matrix
L = g.laplacian();

% Sparse Laplacian with weights from the "R" atribute of the nodes
LR = g.laplacian('field', "R");

% Algebraic connectivity
f = g.fiedler();
c = g.isConnected();