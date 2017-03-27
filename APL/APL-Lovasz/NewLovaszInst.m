%%%% the code is used to generate a random connected graph

% For each node we need at least one edge.

% Start with one node. In each iteration, create a new node and a new edge. 
% The edge is to connect the new node with a random node from the previous node set.

% After all nodes are created, create random edges until the total number of edges (n) is fulfilled. 
% Make sure not to create double edges (for this we use the adjancy matrix).

% Random graph is done in O(n).

function [A,n,ExA] = NewLovaszInst(m, arc)

%%% output: A-the adjacy matrix
%           ExA-the coefficent matrix used in the Lovasz capacity problem
%           n-the actual number of arcs or the dimension of decision
%           variable

%arc = 100; % disgned number of arcs
%m = 10; % number of nodes
n = 0; % actual number of arcs

A = sparse(m,m); %% the adjancy matrix
ExA = sparse(m^2,1);

for i = 2:m,
    %randomly choose a node from the previous node set
    r = ceil((i-1)*rand());
    A(i,r) = 1;
    A(r,i) = 1;
    n = n+1;
    temp = sparse(m^2,1);
    temp((r-1)*m + i) = 1;
    temp((i-1) *m + r) = 1;
    ExA = [ExA, temp];
end

for i = 1:arc,
    r1 = ceil(m*rand());
    
    prop = rand();
    if r1 > 1 && prop < r1 / m,
        r2 = ceil((r1-1) * rand());
    else
        r2 = r1 + ceil((m-r1)*rand());
    end
    if A(r1,r2) == 0
        A(r1,r2) = 1;
        A(r2,r1) = 1;
        temp = sparse(m^2,1);
        temp((r1-1)*m + r2) = 1;
        temp((r2-1) *m + r1) = 1;
        ExA = [ExA, temp];
        n = n+1;
    end
end

disp(sprintf('actual number of edges %d',n));

% the original obejctive function is given by
% \lambda_max(d + \sum_{i<j} x_{ij} (E_{ij} + E_{ij}))
% and d = 1 - A
% we formulate the extended A as follows
d = 1-A;
ExA(:,1) = reshape(d, [],1);


