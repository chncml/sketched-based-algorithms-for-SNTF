function test_different_resample_fixed
clc; clear; close all
format short e
I = 400;J = 20;R = J;

maxIter = 1e10;maxTime = 500;tol = 1e-5;
computeobj = false;

 % generating the syntheic data
Gtrue = rand(I, J);Ytensor = ktensor({Gtrue, Gtrue, Gtrue});
Y = double(tensor(Ytensor));

% N = randn(I, I, I);N = symmetrize(tensor(N));N = double(N);sN = norm(N(:));
% sY = norm(Y(:));ratio = sY/(sN * sqrt(10^(SNR/10)));
% Y = max(Y + N * ratio, 0);
 
% Monte Carlo tests


MCT = 1; NumofMCT = 10; 
while MCT <= NumofMCT
    
    MCT
    
    % initialization
    G0 = rand(I, R) + 1e-5;
%     G0 = ones(I, R);
    
    % performing the proposed multiplicative algorithms
    tic;[G1, f1, t1, fit1] = rand_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G2, f2, t2, fit2] = rand_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G3, f3, t3, fit3] = randkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G4, f4, t4, fit4] = randkr_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G5, f5, t5, fit5] = uniform_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, computeobj);toc
    tic;[G6, f6, t6, fit6] = uniform_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, computeobj);toc
    tic;[G7, f7, t7, fit7] = uniformkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, computeobj);toc
    tic;[G8, f8, t8, fit8] = uniformkr_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, computeobj);toc
    
    T1(MCT,1) = t1(end);
    T2(MCT,1) = t2(end);
    T3(MCT,1) = t3(end);
    T4(MCT,1) = t4(end);
    T5(MCT,1) = t5(end);
    T6(MCT,1) = t6(end);
    T7(MCT,1) = t7(end);
    T8(MCT,1) = t8(end);


    Fit1(MCT,1) = fit1(end);
    Fit2(MCT,1) = fit2(end);
    Fit3(MCT,1) = fit3(end);
    Fit4(MCT,1) = fit4(end);
    Fit5(MCT,1) = fit5(end);
    Fit6(MCT,1) = fit6(end);
    Fit7(MCT,1) = fit7(end);
    Fit8(MCT,1) = fit8(end);

    MCT = MCT+1;

end
        
disp('The results are shown as follows:')

meanFIT = [mean(Fit1),mean(Fit2),mean(Fit3),mean(Fit4),mean(Fit5),mean(Fit6),mean(Fit7),mean(Fit8)]
stdFIT = [std(Fit1),std(Fit2),std(Fit3),std(Fit4),std(Fit5),std(Fit6),std(Fit7),std(Fit8)]

meanT = [mean(T1),mean(T2),mean(T3),mean(T4),mean(T5),mean(T6),mean(T7),mean(T8)]
stdT = [std(T1),std(T2),std(T3),std(T4),std(T5),std(T6),std(T7),std(T8)]


I = 400;J = 20;R = J;SNR = 10; 

maxIter = 1e10;maxTime = 500;tol = 1e-5;
computeobj = false;

 % generating the syntheic data
Gtrue = rand(I, J);Ytensor = ktensor({Gtrue, Gtrue, Gtrue});
Y = double(tensor(Ytensor));

N = randn(I, I, I);N = symmetrize(tensor(N));N = double(N);sN = norm(N(:));
sY = norm(Y(:));ratio = sY/(sN * sqrt(10^(SNR/10)));
Y = max(Y + N * ratio, 0);
 
% Monte Carlo tests


MCT = 1; NumofMCT = 10; 
while MCT <= NumofMCT
    
    MCT
    
    % initialization
    G0 = rand(I, R) + 1e-5;
%     G0 = ones(I, R);
    
    % performing the proposed multiplicative algorithms
    tic;[G1, f1, t1, fit1] = rand_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G2, f2, t2, fit2] = rand_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G3, f3, t3, fit3] = randkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G4, f4, t4, fit4] = randkr_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc
    tic;[G5, f5, t5, fit5] = uniform_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, computeobj);toc
    tic;[G6, f6, t6, fit6] = uniform_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 0.25, computeobj);toc
    tic;[G7, f7, t7, fit7] = uniformkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, computeobj);toc
    tic;[G8, f8, t8, fit8] = uniformkr_Parallel_Multi_SNTF_resample(Y, G0, (1/2), maxIter, maxTime, tol, 0.5, computeobj);toc
    
    T1(MCT,1) = t1(end);
    T2(MCT,1) = t2(end);
    T3(MCT,1) = t3(end);
    T4(MCT,1) = t4(end);
    T5(MCT,1) = t5(end);
    T6(MCT,1) = t6(end);
    T7(MCT,1) = t7(end);
    T8(MCT,1) = t8(end);


    Fit1(MCT,1) = fit1(end);
    Fit2(MCT,1) = fit2(end);
    Fit3(MCT,1) = fit3(end);
    Fit4(MCT,1) = fit4(end);
    Fit5(MCT,1) = fit5(end);
    Fit6(MCT,1) = fit6(end);
    Fit7(MCT,1) = fit7(end);
    Fit8(MCT,1) = fit8(end);

    MCT = MCT+1;

end
        
disp('The results are shown as follows:')

meanFIT = [mean(Fit1),mean(Fit2),mean(Fit3),mean(Fit4),mean(Fit5),mean(Fit6),mean(Fit7),mean(Fit8)]
stdFIT = [std(Fit1),std(Fit2),std(Fit3),std(Fit4),std(Fit5),std(Fit6),std(Fit7),std(Fit8)]

meanT = [mean(T1),mean(T2),mean(T3),mean(T4),mean(T5),mean(T6),mean(T7),mean(T8)]
stdT = [std(T1),std(T2),std(T3),std(T4),std(T5),std(T6),std(T7),std(T8)]
end
