%ENME 547 Final Project F23
%Cyril Justin UCID:30070873

clc;

%Preprocessing

%Discretizing

Lx = 1; %length of square in x direction.
Ly = 1; %Length of square in y direction.
n = input("Enter number of elements: ");
xNodes = linspace(0, Lx, n+1);
yNodes = linspace(0,Ly,n+1);
[X,Y] = meshgrid(xNodes,yNodes);

% Reshape the coordinates to a column vector
nodal_coord = [X(:), Y(:)];

%Define number of elements and number of nodes
nel=n^2;
nen=(n+1)^2;

%Initialize arrays ID, IEN and LM

ID = zeros([nen,1]);
IEN = zeros(4,nel);
LM = zeros(4,nel);

% Assign values to ID array
count = 1;
for i = 1:nen
    if nodal_coord(i, 1) == 0 || nodal_coord(i, 2) == 0
        ID(i) = 0;  % Nodes on x or y axis
    else
        ID(i) = count;  % Other nodes
        count = count + 1;
    end
end
disp('ID Array:');
disp(ID);

% Populate IEN array
k = 1;
for j = 1:n
    for i = 1:n
        node1 = i + (j - 1) * (n + 1);
        node2 = i + 1 + (j - 1) * (n + 1);
        node3 = i + 1 + j * (n + 1);
        node4 = i + j * (n + 1);

        IEN(:, k) = [node1; node2; node3; node4];
        k = k + 1;
    end
end
disp('IEN Array:');
disp(IEN);

% Populate LM array based on ID array
for i = 1:4
    for j = 1:nel
        LM(i, j) = ID(IEN(i, j));
    end
end

% Display the LM array
disp('LM Array:');
disp(LM);
