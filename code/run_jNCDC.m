data=load('..\data\synfix_3.mat');
data=cat(2,data.data);
Label=load('..\data\synfix_3_label.mat');
label=cat(2,Label.label);
[n,~,T]=size(data);
pNmis = zeros(10, T-1);
Nmis = zeros(10, T);
Accs = zeros(10, T);
Fs = zeros(10, T);
RIs = zeros(10, T);
pmeans_nmi = zeros(1, T-1);
means_nmi = zeros(1, T);
means_acc = zeros(1, T);
means_f = zeros(1, T);
means_ri = zeros(1, T);

dim = 80;
alpha = 0.1;
beta = 0.001;
seta1 = 0.8;
seta2 = 0.8;
iter = 50;
comm = 4;
for i = 1:10
    [Z, F, err] = jNCDC(data, iter, dim, alpha, beta, seta1, seta2);
    Y = zeros(T, n);
    C = zeros(n, n);
    for t = 1:T
        if t == 1 || t == T
            C=F(:,:,t)*F(:,:,t)';
        else
           C = (Z(:,:,t)+Z(:,:,t)')/2;
         end
        Y(t,:) = kmeans(C, comm, 'Replicates',10);
        [result, res] = ClusteringMeasure_new(label(t,:), Y(t,:));
        Accs(i,t) = result(1, 1);
        Nmis(i,t) = result(1,2);
        Fs(i,t) = result(1, 3);
        RIs(i,t) = result(1, 4);
        if t ~= 1
            pNmis(i, t-1) = NMI(Y(t-1,:), Y(t,:));
        end
    end
       
end

for t=1:T
   means_nmi(1,t)=mean(Nmis(:,t));
   means_acc(1,t)=mean(Accs(:,t));
   means_f(1,t)=mean(Fs(:,t));
   means_ri(1,t)=mean(RIs(:,t));
end
ave_nmi = mean(means_nmi);
ave_acc = mean(means_acc);
ave_f = mean(means_f);
ave_ri = mean(means_ri);


disp(["jNCDC synfix_AVE_NMI: ", num2str(ave_nmi)]);
disp(["ave_nmi", num2str(ave_nmi)]);
disp(["ave_acc", num2str(ave_acc)]);
disp(["ave_f", num2str(ave_f)]);
disp(["ave_ri", num2str(ave_ri)]);






