function [Z,F,err] = jNCDC(W,Iter,k,alpha,beta,seta1,seta2)

if nargin>8
	error('parameter is too much,the max number of parameter is 5');
end
switch nargin
	case 1
		Iter=100;
		k=80;
		alpha=0.1;
		beta=1;
	case 2
		k=80;
		alpha=0.1;
		beta=1;
	case 3 
		alpha=0.1;
		beta=1;
	case 4
		beta=1;
end

[n,~,T]=size(W);

%Initialization M,L,D
M=zeros(n,n,T);
D=zeros(n,n,T);

%Caculation M,L,D
for i=1:T
    M(:,:,i)=PMI(W(:,:,i),2);  
end

%L=LMatrix(W);

for i=1:T
    D(:,:,i)=diag(sum(W(:,:,i),2));
end 
% L = D-W;
%Initialize err
err=zeros(Iter,1,T);

%Initialize dindex,sindex;
[dindex,sindex]=Divide(W,0.05);

% 初始化B,F
B=zeros(n,k);
F=zeros(n,k,T);
Z=zeros(n,n,T);
for o=1:T
    pM=M(:,:,o);
    [b,c,f]=svds(pM,k);
    F(:,:,o)=abs((c^0.5)*f')';
    F(:,:,o)=mapminmax(F(:,:,o),0,1);
    if o==1
        B=abs(b*(c^0.5));
        B=mapminmax(B,0,1);
    elseif o == T
        B_10 = abs(b*(c^0.5));
        B_10=mapminmax(B_10,0,1);
    else
        [fu,fs,fv]=svds(F(:,:,o),k,'largest','SubspaceDimension',200);
        Z(:,:,o)=abs((fu*(fs^0.5))*(fu*(fs^0.5))');
        Z(:,:,o)=mapminmax(Z(:,:,o),0,1);
        clear fu;
        clear fs;
        clear fv;
    end
	clear b;
    clear c;
end
A=B;
G=Z;
E=F;

for i=1:T
    % Adam
    averageGrad = [];
    averageSqGrad = [];
    learnRate = 0.05;
%     gradDecay = 0.75;
%     sqGradDecay = 0.95;
%     learnRate = 0.001;
    gradDecay = 0.9;
    sqGradDecay = 0.999;
    iteration=3;
    if i==1
        L1 = D(:,:,i) - W(:,:,i);
        for o=1:Iter  
            B(:,:)=A(:,:).*((M(:,:,i)*E(:,:,i))./(A(:,:)*E(:,:,i)'*E(:,:,i)+eps));
            F(:,:,i)=E(:,:,i).*((M(:,:,i)'*A(:,:)+W(:,:,i)*E(:,:,i))./(E(:,:,i)*A(:,:)'*A(:,:)+D(:,:,i)*E(:,:,i)+eps));
            B(:,:)=mapminmax(B(:,:),0,1);
            F(:,:,i)=mapminmax(F(:,:,i),0,1);
            err(o,1,i)=norm(M(:,:,1)-B*F(:,:,1)','fro')^2 + trace(F(:,:,i)'*L1*F(:,:,i));
            order=log10(err(o,1,i));
            order = 1* 10^(-order);
            if o~=1
                if abs((err(o-1,1,i)-err(o,1,i))/err(o,1,i)) <= order*0.1
                    break
                else
                    A(:,:)=B(:,:);
                    E(:,:,1)=F(:,:,1);
                end
            else
                A(:,:)=B(:,:);
                E(:,:,1)=F(:,:,1);
            end
        end
        [fu,fs,fv]=svds(F(:,:,1));
        Z(:,:,i)=abs((fu*(fs^0.5))*((fu*(fs^0.5)))');
        Z(:,:,i)=mapminmax(Z(:,:,i),0,1);
        G(:,:,i) = Z(:,:,i);
        clear fu;
        clear fs;
        clear fv;
        clear L1;
    elseif i == T
        B = B_10;
        clear B_10;
        A = B;
        L1 = D(:,:,i) - W(:,:,i);
        for o=1:Iter
            B(:,:)=A(:,:).*((M(:,:,i)*E(:,:,i))./(A(:,:)*E(:,:,i)'*E(:,:,i)+eps));
            F(:,:,i)=E(:,:,i).*((M(:,:,i)'*A(:,:)+W(:,:,i)*E(:,:,i))./(E(:,:,i)*A(:,:)'*A(:,:)+D(:,:,i)*E(:,:,i)+eps));
            B(:,:)=mapminmax(B(:,:),0,1);
            F(:,:,i)=mapminmax(F(:,:,i),0,1);
            err(o,1,i)=norm(M(:,:,i)-B*F(:,:,i)','fro')^2+ trace(F(:,:,i)'*L1*F(:,:,i));
            order=log10(err(o,1,i));
            order = 1* 10^(-order);
            if o~=1
                if abs((err(o-1,1,i)-err(o,1,i))/err(o,1,i)) <= order*0.1
                    break
                else
                    A(:,:)=B(:,:);
                    E(:,:,i)=F(:,:,i);
                end
            else
                A(:,:)=B(:,:);
                E(:,:,i)=F(:,:,i);
            end
        end
        clear L1;
    else
        
        SMF=(1/3)*(M(:,:,i-1)*E(:,:,i-1)+M(:,:,i)*E(:,:,i)+M(:,:,i+1)*E(:,:,i+1));
        [U,S,V]=svds(SMF);
        clear SMF;
        B=abs(U*V');
        B(find(isnan(B)==1)) = 0;
        B = mapminmax(B,0,1);
        A=B;
        nowDindex=dindex(:,:,i-1);
        nowSindex=sindex(:,:,i-1);
        if i~=2
            preDindex=dindex(:,:,i-2);
            preSindex=sindex(:,:,i-2);
        end
        if i~=T
            NextDindex=dindex(:,:,i);
            NextSindex=sindex(:,:,i);
        end
        for o=1:Iter
            grad=zeros(n,n);              
            B = A.*((M(:,:,i-1)*E(:,:,i-1) + M(:,:,i)*E(:,:,i) + M(:,:,i+1)*E(:,:,i+1))./(A*(E(:,:,i-1)'*E(:,:,i-1)+E(:,:,i)'*E(:,:,i)+E(:,:,i+1)'*E(:,:,i+1))+eps));
            F(:,:,i)=E(:,:,i).*((M(:,:,i)'*A+alpha*G(:,:,i)*E(:,:,i-1)+alpha*G(:,:,i+1)'*E(:,:,i+1)+W(:,:,i)*E(:,:,i))./(E(:,:,i)*(A')*A+D(:,:,i)*E(:,:,i)+alpha*E(:,:,i)+alpha*G(:,:,i+1)*G(:,:,i+1)'*E(:,:,i) +eps));
            TF=F(:,:,i);
            TF(find(isnan(TF)==1)) = 0;
            %F(:,:,i) = TF;
            B = mapminmax(B, 0, 1);
            F(:,:,i) = mapminmax(TF, 0, 1);
            %E(:,:,i) = F(:,:,i);

            clear TF;
            [fu,fs,fv]=svds(F(:,:,i),k,'largest','SubspaceDimension',100); % 100
            sz = abs((fu*fs^0.5)*(fu*fs^0.5)');
            sz(find(isnan(sz)==1))=0;
            Z(:,:,i) = mapminmax(sz,0,1);
            clear sz;
            clear fu;
            clear fs;
            clear fv;
            currentDNM = findNeighbor_3(nowDindex, G(:,:,i), seta1);
            currentDSM = findNeighbor_3(nowSindex, Z(:,:,i-1),seta2);
            if i ~= 2
               PastDNM = findNeighbor_3(preDindex, Z(:,:,i-1), seta1);
               PastSNM = findNeighbor_3(preSindex, Z(:,:,i-2), seta2);
            end
            NextDNM = findNeighbor_3(NextDindex, Z(:,:,i+1), seta1);
            NextSNM = findNeighbor_3(NextSindex, G(:,:,i),seta2);
             FF=E(:,:,i)*E(:,:,i-1)';     
             ZF=G(:,:,i-1)*E(:,:,i-1)*E(:,:,i-1)';
             DGrad = zeros(n,n);
             SGrad = zeros(n,n);
             DGrad = CaculationGrad(currentDNM, G(:,:,i), FF, ZF, alpha, beta);
             SGrad = CaculationGrad(currentDSM, G(:,:,i), FF, ZF, alpha, beta);
             grad = DGrad + SGrad;
             clear DGrad;
             clear SGrad;
             % update Z
            currentZ=Z(:,:,i);
            [currentZ,averageGrad,averageSqGrad] = adamupdate(currentZ,grad,averageGrad,averageSqGrad,iteration,learnRate,gradDecay,sqGradDecay);
            currentZ(find(currentZ<0))=0;
            currentZ(find(isnan(currentZ)==1))=0;
            Z(:,:,i) = mapminmax(currentZ, 0, 1);
            clear currentZ;
           
            if i==2
                L1 = D(:,:,i-1) - W(:,:,i-1);
                L2 = D(:,:,i) - W(:,:,i);
                L3 = D(:,:,i+1) - W(:,:,i+1);
                loss1=norm(M(:,:,i-1)-B(:,:,1)*F(:,:,i-1)','fro')^2+norm(M(:,:,i)-B(:,:,1)*F(:,:,i)','fro')^2+norm(M(:,:,i+1)-B(:,:,1)*F(:,:,i+1)','fro')^2;
                loss2=trace(F(:,:,i)'*L2*F(:,:,i));
                loss3=norm(F(:,:,i)-Z(:,:,i)*F(:,:,i-1),'fro')^2;
                loss4_2 = CaculationContrastiveLoss_3(currentDNM, Z(:,:,i)) + CaculationContrastiveLoss_3(currentDSM, Z(:,:,i));
                loss4_3 = CaculationContrastiveLoss_3(NextDNM, Z(:,:,i+1)) + CaculationContrastiveLoss_3(NextSNM, Z(:,:,i+1));
                loss4=loss4_2+loss4_3;
                clear currentDNM;
                clear currentDSM;
                clear NextDNM;
                clear NextSNM;
                clear L1;
                clear L2;
                clear L3;
            elseif i==T         
                L1 = D(:,:,i-1) - W(:,:,i-1);
                L2 = D(:,:,i) - W(:,:,i);
                loss1=norm(M(:,:,i-1)-B(:,:,1)*F(:,:,i-1)','fro')^2+norm(M(:,:,i)-B(:,:,1)*F(:,:,i)','fro')^2;
                loss2=+trace(F(:,:,i)'*L2*F(:,:,i));
                loss3=norm(F(:,:,i)-Z(:,:,i)*F(:,:,i-1),'fro')^2;
                loss4_2 = CaculationContrastiveLoss_3(currentDNM, Z(:,:,i)) + CaculationContrastiveLoss_3(currentDSM, Z(:,:,i));
                loss4_1 = CaculationContrastiveLoss_3(PastDNM, Z(:,:,i-1)) + CaculationContrastiveLoss_3(PastSNM, Z(:,:,i-1));
                loss4=loss4_2+loss4_1;
                clear currentDNM;
                clear currentDSM;
                clear PastDNM;
                clear PastSNM;
                clear L1;
                clear L2;
                % clear L3;
            else
                L1 = D(:,:,i-1) - W(:,:,i-1);
                L2 = D(:,:,i) - W(:,:,i);
                L3 = D(:,:,i+1) - W(:,:,i+1);
                loss1=norm(M(:,:,i-1)-B(:,:,1)*F(:,:,i-1)','fro')^2+norm(M(:,:,i)-B(:,:,1)*F(:,:,i)','fro')^2+norm(M(:,:,i+1)-B(:,:,1)*F(:,:,i+1)','fro')^2;
                loss2=trace(F(:,:,i)'*L2*F(:,:,i));
                loss3=norm(F(:,:,i)-Z(:,:,i)*F(:,:,i-1),'fro')^2;
                loss4_1 = CaculationContrastiveLoss_3(PastDNM, Z(:,:,i-1)) + CaculationContrastiveLoss_3(PastSNM, Z(:,:,i-1));
                loss4_2 = CaculationContrastiveLoss_3(currentDNM, Z(:,:,i)) + CaculationContrastiveLoss_3(currentDSM, Z(:,:,i));
                loss4_3 = CaculationContrastiveLoss_3(NextDNM, Z(:,:,i+1)) + CaculationContrastiveLoss_3(NextSNM, Z(:,:,i+1));
                loss4=loss4_2+loss4_1+loss4_3;
                clear currentDNM;
                clear currentDSM;
                clear PastDNM;
                clear PastSNM;
                clear NextDNM;
                clear NextSNM;
                clear L1;
                clear L2;
                clear L3;
            end
            
            err(o,1,i)=loss1+ loss2+ alpha*loss3+ beta*loss4; 
            if o~=1
                deloss = (err(o-1,1,i)-err(o,1,i))/err(o,1,i);
               if abs(deloss)<1e-7
                   break
               else
                   A=B;
                   E(:,:,i)=F(:,:,i);
                   G(:,:,i)=Z(:,:,i);
               end
            else
                A=B;
                E(:,:,i)=F(:,:,i);
                G(:,:,i)=Z(:,:,i);
            end
         end
    end
 end
end