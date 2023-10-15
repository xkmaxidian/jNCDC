function [loss] = CaculationContrastiveLoss_3(NM, Z)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
%  --Inputs
%    NM:A [N, N] Matrix, store neighbors
%    Z: A [n, n] Matrix, simility matrix;
%  --Outputs
%    loss: A scalar
% Initialize losaa
loss = 0;
diags = diag(Z, 0);                                % 获取矩阵的对角线元素
k0 = sum(exp(Z),2)-exp(diags);                  % 数组版：k0=sum(exp(Z(nd_i,:)),2)-exp(Z(nd_i,nd_i));
NM(logical(eye(size(NM))))=0;                   % 将邻居矩阵的对角线元素置为0，if nd_i~=nd_j，防止计算对角线元素
% 计算每个节点的loss,即1行
lossRows = -log10(sum(exp(NM.*Z),2)./k0);            %-log10(exp(Z(nd_i,nd_j))/k0);
loss = sum(lossRows, 1);
end