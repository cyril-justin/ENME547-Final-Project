%ENME 547 Final Project F23
%Cyril Justin UCID:30070873

clc;clear;

%Preprocessing

%Discretizing

Lx = 1; %length of square in x direction.
Ly = 1; %Length of square in y direction.
global_coord = [0,0;1,0;1,1;0,1];
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

%Defining the basis functions

ksi = input('Please enter the ksi coordinate (between -1 and 1): ');
eta = input('Please enter the eta coordinate (between -1 and 1): ');


Na1 = (1/4)*(1-ksi)*(1-eta)
Na2 = (1/4)*(1-ksi)*(1-eta);
Na3 = (1/4)*(1-ksi)*(1-eta);
Na4 = (1/4)*(1-ksi)*(1-eta);

    dNa1_dxi = @(ksi, eta) -0.25*(1-eta);
    dNa2_dxi = @(ksi, eta) 0.25*(1-eta);
    dNa3_dxi = @(ksi, eta) 0.25*(1+eta);
    dNa4_dxi = @(ksi, eta) -0.25*(1+eta);

    dNa1_deta = @(ksi, eta) -0.25*(1-ksi);
    dNa2_deta = @(ksi, eta) -0.25*(1+ksi);
    dNa3_deta = @(ksi, eta) 0.25*(1+ksi);
    dNa4_deta = @(ksi, eta) 0.25*(1-ksi);

    DDDNa1=dNa1_dxi(ksi,eta)
    DDNa1=diff(Na1,ksi)  

