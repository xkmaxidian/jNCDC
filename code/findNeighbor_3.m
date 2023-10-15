function [NeighborMatrix] = findNeighbor_3(index,Z,seta)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
% --Input
%   --index, 动态节点的索引;
%   --Z is a [n,n] matrix, n is node number, similarity matrix;
%   --seta, 阈值

% --output
%   --NeighborMatrix [n, n] 动态节点的邻居节点矩阵，如NM[i, j]=1,说明j是i的邻居节点，反之不为邻居

% 首先按每一行排序
[n, ~] = size(Z);
knodes = 15;
SMatrix = zeros(n, n);
[SMatrix, node_index] = sort(Z, 2, 'descend');

clear SMatrix;
NeighborMatrix = zeros(n, n);
%  NeighborMatrix(index, node_index(:,1:10)) = 1;
NeighborMatrix(index, node_index(:,1: knodes)) = 1;
end