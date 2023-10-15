function [Grad] = CaculationGrad(NM, Z, FF, ZF, alpha, beta)
%UNTITLED5 此处提供此函数的摘要
%   此处提供详细说明
%  --Inputs
%    NM: A [n, n] matrix, store neigbors;
%    Z : A [n, n] matrix, simility matrix;
%    FF: A [n, n] matrix, F(:,:)*F(:,:)';
%    ZF: A [n, n] matrix, Z(:,:)*F(:,:,i-1)*F(:,:,i-1)';
%    alpha: A scalar, parameter
%    beta: A scalar, parameter

%   --outputs
%    Grad：A [n, n] matrix; Grad matrix;

[n, ~] = size(Z);
Grad = zeros(n, n);
F1 = zeros(n,n);

diags = diag(Z, 0);
k0 = sum(exp(Z),2) - exp(diags);
NM(logical(eye(size(NM)))) = 0;
knum = sum(NM, 2);                 
F1 = -2 * FF + 2 * ZF;
tempNM1 = NM;

tempNM1(find(tempNM1==1))=-1;
tempNM2 = tempNM1 + 1;

F2 = (tempNM1 + knum.*(sum(exp(NM.*Z),2)./k0)) + (knum.*(sum(exp(tempNM2.*Z),2)./k0));
Grad = alpha * F1 + beta *  F2;

end