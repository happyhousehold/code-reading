%%%% the code is used to generate a random connected graph

% For each node we need at least one edge.

% Start with one node. In each iteration, create a new node and a new edge. 
% The edge is to connect the new node with a random node from the previous node set.

% After all nodes are created, create random edges until the total number of edges (n) is fulfilled. 
% Make sure not to create double edges (for this we use the adjancy matrix).

% Random graph is done in O(n).

function [A,narc] = NewLovaszInst(m, arc)


%arc = 100; % disgned number of arcs
%m = 10; % number of nodes
narc = 0; % actual number of arcs

A = sparse(m,m); %% the adjancy matrix

for i = 2:m,
    %randomly choose a node from the previous node set
    r = ceil((i-1)*rand());
    A(i,r) = 1;
    A(r,i) = 1;
    narc = narc+1;
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
        narc = narc+1;
    end
end


