function [dindex,sindex] = Divide(W,top)
% --Input  
%   --W is a [n,n,T] matrix, n is node number, T is the total time
%   --top is a percentage, for select dynamic node


%  对于动态节点，其实只需知道其索引便可


% --Output
%   --dindex is a [nd,1,T-1] matrix, nd is dynamic node number, T is the total
%   time
%   --Wd is a [nd,n,T-1] matrix
%   --sindex is a [ns,1,T-1] matrix, ns is static node number, ns=n-nd
%   --Ws is a [ns,n,T-1] matrix


[n,~,T]=size(W);

%Initialization dindex,sindex,Wd,Ws,nd
nd=round(n*top);         % 强制类型转换，确保结果为整数，round()四舍五入
ns=n-nd;
dindex=zeros(nd,1,T-1);
sindex=zeros(ns,1,T-1);
% Wd=zeros(nd,n,T-1);
% Ws=zeros(ns,n,T-1);
deltaW=zeros(n,n,T-1);
%Caculation deltaW for to find dynamic node
for i=2:T
    % deltaW(:,:,i-1)=W(:,:,i)-W(:,:,i-1); 
    deltaW=W(:,:,i)-W(:,:,i-1); 
    deltaWSum=sum(abs(deltaW),2);          % [n,1,T-1],对deltaW求行和，以获取节点的变化情况
    clear deltaW;
    [deltaWSum,index]=sort(deltaWSum,1,"descend");  %降序排序
    for j=1:nd
        dindex(j,1,i-1)=index(j,1);
        % dindex(j,1,i-1)=index(j,1,i-1);   /2022/10/12
        %midindex=index(j,1,i-1);
        %Wd(j,:,i-1)=W(midindex,:,i-1);
    end
    % 对静态节点的选取做处理，选取10%作为静态节点
    ns = n-nd;
    for k=1:ns
        sindex(k,1,i-1)=index(nd+k,1);
        %sindex(k,1,i-1)=index(nd+k,1,i-1); /2022/10/12
        %midindex=index(nd+k,1,i-1);
        %Ws(k,:,i-1)=W(midindex,:,i-1);
    end
    %nodenum(i-1,1)=deltaWSum(1,1,i-1);
end
end