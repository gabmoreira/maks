%%
clc;
close all;
clear all;

addpath('geo');
addpath('math');

g = Posegraph();

% Specify the attributes of each node (type/shape agnostic)
g.configNodes("R", "t");

% Specify the attributes of each edge (type agnostic)
g.configEdges("R", "t");

% Add nodes with attributes specified earlier
g.addNode(1, 'R', rand(3,3), 't', rand(2,1))
g.addNode(2, 'R', rand(3,3), 't', rand(2,1))

% Add edges with attributes specified earlier between existing nodes
g.addEdge([1;2], 'R', rotrnd(1,0.5), 't', rand(3,1));
g.addEdge([2;1], 'R', g(1,2).R', 't', rand(3,1));

%%

addpath('io');
dataset = 'parking-garage';

[ij, R, t] = readG2O("../data/" + dataset + "_SO3.g2o");


%% Accessing nodes and edges structs

% Node with id 1
g(1)

% Edge from id 1 to id 2
g(1,2)

% Error
g(4,10)

%% Accessing nodes and edges attributes

g(1).R
g(1,2).R

% Error
g(1).T

%% Apply a function to a certain attribute in all the nodes 

g.nodeOp(@projectToSO3, 'R');

%% Graph structure

% As vector
g.degree()

% As sparse diagonal matrix
g.degree('asmatrix', true)

% Sparse adjacency matrix
g.adjacency()

% Sparse Laplacian matrix
g.laplacian()

% Algebraic connectivity
g.fiedler()
g.isConnected()