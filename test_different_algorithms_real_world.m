function test_different_algorithms_real_world
% Experiment on real-world data

clc; clear; close all 

dataset = 'UMIST';
K = 20;

switch dataset
    case {'COIL'}
          load COIL20.mat 
    case {'ORL'}
          load ORL_32x32.mat
    case {'UMIST'}
          load umist.mat 
          fea = fea';
    case {'YALE'}
          load YaleB.mat
          fea = fea';
    otherwise        
          error('the dataset does not exist!');
end    


fulldata = fea';
numperclass = 10;
maxIter = 1e10;
maxTime = 500;
tol = 1e-5;
computeobj = false;
tau = 0.1;
    
%% Monte Carlo tests (MCT)

MCT = 1; NumofMCT = 10;  
while MCT <= NumofMCT
      
    % generate data subset        
    [datasub,subN,labelsub,Trans_datasub,clist] = loadsub(fulldata,gnd,K,numperclass); 
       
    X = datasub;  
    label = labelsub;
           
    [U,~,latent] = pca(X');
    X = U(:,1:K)'*X;
       
   %% construct affinity tensor 
    
    [nfea,ndata] = size(X);

    X2 = sum(X.^2);
    D = repmat(X2,ndata,1)+repmat(X2',1,ndata)-2*X'*X;
    D = real(sqrt(D));
      
    % construct 3-way affinity tensor

    V = zeros([ndata,ndata,ndata]);
    for i = 1:ndata
        i
        for j = i:ndata  
            for k = j:ndata  
                temp = D(i,j)+D(i,k)+D(j,k);
                V(i,j,k) = temp;
                V(i,k,j) = temp;
                V(j,i,k) = temp;
                V(j,k,i) = temp;
                V(k,i,j) = temp;
                V(k,j,i) = temp;
            end
        end
    end
    
    V2 = V.^2;
    temp = sort(V2(:),1,'ascend');
    delta = temp(floor(tau*(ndata^3)))/3;
    Y = exp(-1*V2/delta);
    Y = HyperStochasticTensor(Y);
    Y = max(Y,1E-12);
           
    size(Y)
    %% perform SNTF
    
    for init_iter = 1:5
        
        disp([MCT,init_iter])
        
        G0 = rand(ndata,K)+1e-5; 
        
        tic;[G1, f1, t1, fit1] = Parallel_Multi_SNTF(Y, G0, (1/6), maxIter, maxTime, tol, computeobj);toc;
        tic;[G2, f2, t2, fit2] = Parallel_Multi_SNTF(Y, G0, (1/5), maxIter, maxTime, tol, computeobj);toc;
        tic;[G3, f3, t3, fit3] = rand_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
        tic;[G4, f4, t4, fit4] = randkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
        tic;[G5, f5, t5, fit5] = uniform_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, computeobj);toc
        tic;[G6, f6, t6, fit6] = uniformkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, computeobj);toc
        tic;[G7, f7, t7, fit7] = uniformrand_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, 20, computeobj);toc
        tic;[G8, f8, t8, fit8] = uniformrand_Parallel_Multi_SNTF1(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, 20, computeobj);toc
        tic;[G9, f9, t9, fit9] = uniformrand_Parallel_Multi_SNTF2(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, 20, computeobj);toc
        tic;[G10, f10, t10, fit10] = uniformrand_Parallel_Multi_SNTF3(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, 20, computeobj);toc
        tic;[G_alpha,f_alpha,t_alpha,fit_alpha] = alpha_Paralle_Multi_SNTF(Y,G0,(1/2),maxIter,maxTime,tol);toc
        tic;[G_beta,f_beta,t_beta,fit_beta] = beta_Paralle_Multi_SNTF(Y,G0,(1/2),maxIter,maxTime,tol);toc
        
        [AC_SNTF1(init_iter,MCT),MIhat_SNTF1(init_iter,MCT)] = AC_MIhat(G1,label,K);
        [AC_SNTF2(init_iter,MCT),MIhat_SNTF2(init_iter,MCT)] = AC_MIhat(G2,label,K);
        [AC_SNTF3(init_iter,MCT),MIhat_SNTF3(init_iter,MCT)] = AC_MIhat(G3,label,K);
        [AC_SNTF4(init_iter,MCT),MIhat_SNTF4(init_iter,MCT)] = AC_MIhat(G4,label,K);
        [AC_SNTF5(init_iter,MCT),MIhat_SNTF5(init_iter,MCT)] = AC_MIhat(G5,label,K);
        [AC_SNTF6(init_iter,MCT),MIhat_SNTF6(init_iter,MCT)] = AC_MIhat(G6,label,K);
        [AC_SNTF7(init_iter,MCT),MIhat_SNTF7(init_iter,MCT)] = AC_MIhat(G7,label,K);
        [AC_SNTF8(init_iter,MCT),MIhat_SNTF8(init_iter,MCT)] = AC_MIhat(G8,label,K);
        [AC_SNTF9(init_iter,MCT),MIhat_SNTF9(init_iter,MCT)] = AC_MIhat(G9,label,K);
        [AC_SNTF10(init_iter,MCT),MIhat_SNTF10(init_iter,MCT)] = AC_MIhat(G10,label,K);
        [AC_SNTF_alpha(init_iter,MCT),MIhat_SNTF_alpha(init_iter,MCT)] = AC_MIhat(G_alpha,label,K);
        [AC_SNTF_beta(init_iter,MCT),MIhat_SNTF_beta(init_iter,MCT)] = AC_MIhat(G_beta,label,K);       
      
        T_SNTF1(init_iter,MCT) = t1(end);
        T_SNTF2(init_iter,MCT) = t2(end);
        T_SNTF3(init_iter,MCT) = t3(end);
        T_SNTF4(init_iter,MCT) = t4(end);
        T_SNTF5(init_iter,MCT) = t5(end);
        T_SNTF6(init_iter,MCT) = t6(end);
        T_SNTF7(init_iter,MCT) = t7(end);
        T_SNTF8(init_iter,MCT) = t8(end);
        T_SNTF9(init_iter,MCT) = t9(end);
        T_SNTF10(init_iter,MCT) = t10(end);
        T_SNTF_alpha(init_iter,MCT) = t_alpha(end);
        T_SNTF_beta(init_iter,MCT) = t_beta(end);
                         
    end
            
    MCT = MCT+1;
    
end

[bestAC_SNTF1, index_SNTF1] = max(AC_SNTF1);
bestMIhat_SNTF1 = MIhat_SNTF1((0:NumofMCT-1)*5+index_SNTF1);
bestT_SNTF1 = T_SNTF1((0:NumofMCT-1)*5+index_SNTF1);

[bestAC_SNTF2, index_SNTF2] = max(AC_SNTF2);
bestMIhat_SNTF2 = MIhat_SNTF2((0:NumofMCT-1)*5+index_SNTF2);
bestT_SNTF2 = T_SNTF2((0:NumofMCT-1)*5+index_SNTF2);

[bestAC_SNTF3, index_SNTF3] = max(AC_SNTF3);
bestMIhat_SNTF3 = MIhat_SNTF3((0:NumofMCT-1)*5+index_SNTF3);
bestT_SNTF3 = T_SNTF3((0:NumofMCT-1)*5+index_SNTF3);

[bestAC_SNTF4, index_SNTF4] = max(AC_SNTF4);
bestMIhat_SNTF4 = MIhat_SNTF4((0:NumofMCT-1)*5+index_SNTF4);
bestT_SNTF4 = T_SNTF4((0:NumofMCT-1)*5+index_SNTF4);

[bestAC_SNTF5, index_SNTF5] = max(AC_SNTF5);
bestMIhat_SNTF5 = MIhat_SNTF5((0:NumofMCT-1)*5+index_SNTF5);
bestT_SNTF5 = T_SNTF5((0:NumofMCT-1)*5+index_SNTF5);

[bestAC_SNTF6, index_SNTF6] = max(AC_SNTF6);
bestMIhat_SNTF6 = MIhat_SNTF6((0:NumofMCT-1)*5+index_SNTF6);
bestT_SNTF6 = T_SNTF6((0:NumofMCT-1)*5+index_SNTF6);

[bestAC_SNTF7, index_SNTF7] = max(AC_SNTF7);
bestMIhat_SNTF7 = MIhat_SNTF7((0:NumofMCT-1)*5+index_SNTF7);
bestT_SNTF7 = T_SNTF7((0:NumofMCT-1)*5+index_SNTF7);

[bestAC_SNTF8, index_SNTF8] = max(AC_SNTF8);
bestMIhat_SNTF8 = MIhat_SNTF8((0:NumofMCT-1)*5+index_SNTF8);
bestT_SNTF8 = T_SNTF8((0:NumofMCT-1)*5+index_SNTF8);

[bestAC_SNTF9, index_SNTF9] = max(AC_SNTF9);
bestMIhat_SNTF9 = MIhat_SNTF9((0:NumofMCT-1)*5+index_SNTF9);
bestT_SNTF9 = T_SNTF9((0:NumofMCT-1)*5+index_SNTF9);

[bestAC_SNTF10, index_SNTF10] = max(AC_SNTF10);
bestMIhat_SNTF10 = MIhat_SNTF10((0:NumofMCT-1)*5+index_SNTF10);
bestT_SNTF10 = T_SNTF10((0:NumofMCT-1)*5+index_SNTF10);

[bestAC_SNTF_alpha, index_SNTF_alpha] = max(AC_SNTF_alpha);
bestMIhat_SNTF_alpha = MIhat_SNTF_alpha((0:NumofMCT-1)*5+index_SNTF_alpha);
bestT_SNTF_alpha = T_SNTF_alpha((0:NumofMCT-1)*5+index_SNTF_alpha);

[bestAC_SNTF_beta, index_SNTF_beta] = max(AC_SNTF_beta);
bestMIhat_SNTF_beta = MIhat_SNTF_beta((0:NumofMCT-1)*5+index_SNTF_beta);
bestT_SNTF_beta = T_SNTF_beta((0:NumofMCT-1)*5+index_SNTF_beta);

%% show the results

AC_SNTF = [bestAC_SNTF1',bestAC_SNTF2',bestAC_SNTF_alpha',bestAC_SNTF_beta'...
    ,bestAC_SNTF3',bestAC_SNTF4',bestAC_SNTF5',bestAC_SNTF6'...
    ,bestAC_SNTF7',bestAC_SNTF8',bestAC_SNTF9',bestAC_SNTF10'];
[mean(AC_SNTF);std(AC_SNTF)]

MIhat_SNTF = [bestMIhat_SNTF1',bestMIhat_SNTF2',bestMIhat_SNTF_alpha',bestMIhat_SNTF_beta'...
    ,bestMIhat_SNTF3',bestMIhat_SNTF4',bestMIhat_SNTF5',bestMIhat_SNTF6'...
    ,bestMIhat_SNTF7',bestMIhat_SNTF8',bestMIhat_SNTF9',bestMIhat_SNTF10'];
[mean(MIhat_SNTF);std(MIhat_SNTF)]

T_SNTF = [bestT_SNTF1',bestT_SNTF2',bestT_SNTF_alpha',bestT_SNTF_beta'...
    ,bestT_SNTF3',bestT_SNTF4',bestT_SNTF5',bestT_SNTF6'...
    ,bestT_SNTF7',bestT_SNTF8',bestT_SNTF9',bestT_SNTF10'];
[mean(T_SNTF);std(T_SNTF)]
end