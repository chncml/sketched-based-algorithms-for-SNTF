function test_different_algorithms_comparison
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
    
    T1(MCT,1) = t1(end);
    T2(MCT,1) = t2(end);
    T3(MCT,1) = t3(end);
    T4(MCT,1) = t4(end);
    T5(MCT,1) = t5(end);
    T6(MCT,1) = t6(end);
    T7(MCT,1) = t7(end);
    T8(MCT,1) = t8(end);
    T9(MCT,1) = t9(end);
    T10(MCT,1) = t10(end);
    T_alpha(MCT,1) = t_alpha(end);
    T_beta(MCT,1) = t_beta(end);   


    Fit1(MCT,1) = fit1(end);
    Fit2(MCT,1) = fit2(end);
    Fit3(MCT,1) = fit3(end);
    Fit4(MCT,1) = fit4(end);
    Fit5(MCT,1) = fit5(end);
    Fit6(MCT,1) = fit6(end);
    Fit7(MCT,1) = fit7(end);
    Fit8(MCT,1) = fit8(end);
    Fit9(MCT,1) = fit9(end);
    Fit10(MCT,1) = fit10(end);
    Fit_alpha(MCT,1) = fit_alpha(end);
    Fit_beta(MCT,1) = fit_beta(end);

    MCT = MCT+1;

end
        
disp('The results are shown as follows:')

meanFIT = [mean(Fit1),mean(Fit2),mean(Fit_alpha),mean(Fit_beta),mean(Fit3),mean(Fit4),mean(Fit5),mean(Fit6),mean(Fit7),mean(Fit8),mean(Fit9),mean(Fit10)]
stdFIT = [std(Fit1),std(Fit2),std(Fit_alpha),std(Fit_beta),std(Fit3),std(Fit4),std(Fit5),std(Fit6),std(Fit7),std(Fit8),std(Fit9),std(Fit10)]

meanT = [mean(T1),mean(T2),mean(T_alpha),mean(T_beta),mean(T3),mean(T4),mean(T5),mean(T6),mean(T7),mean(T8),mean(T9),mean(T10)]
stdT = [std(T1),std(T2),std(T_alpha),std(T_beta),std(T3),std(T4),std(T5),std(T6),std(T7),std(T8),std(T9),std(T10)]


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
    
    T1(MCT,1) = t1(end);
    T2(MCT,1) = t2(end);
    T3(MCT,1) = t3(end);
    T4(MCT,1) = t4(end);
    T5(MCT,1) = t5(end);
    T6(MCT,1) = t6(end);
    T7(MCT,1) = t7(end);
    T8(MCT,1) = t8(end);
    T9(MCT,1) = t9(end);
    T10(MCT,1) = t10(end);
    T_alpha(MCT,1) = t_alpha(end);
    T_beta(MCT,1) = t_beta(end);   


    Fit1(MCT,1) = fit1(end);
    Fit2(MCT,1) = fit2(end);
    Fit3(MCT,1) = fit3(end);
    Fit4(MCT,1) = fit4(end);
    Fit5(MCT,1) = fit5(end);
    Fit6(MCT,1) = fit6(end);
    Fit7(MCT,1) = fit7(end);
    Fit8(MCT,1) = fit8(end);
    Fit9(MCT,1) = fit9(end);
    Fit10(MCT,1) = fit10(end);
    Fit_alpha(MCT,1) = fit_alpha(end);
    Fit_beta(MCT,1) = fit_beta(end);

    MCT = MCT+1;

end
        
disp('The results are shown as follows:')

meanFIT = [mean(Fit1),mean(Fit2),mean(Fit_alpha),mean(Fit_beta),mean(Fit3),mean(Fit4),mean(Fit5),mean(Fit6),mean(Fit7),mean(Fit8),mean(Fit9),mean(Fit10)]
stdFIT = [std(Fit1),std(Fit2),std(Fit_alpha),std(Fit_beta),std(Fit3),std(Fit4),std(Fit5),std(Fit6),std(Fit7),std(Fit8),std(Fit9),std(Fit10)]

meanT = [mean(T1),mean(T2),mean(T_alpha),mean(T_beta),mean(T3),mean(T4),mean(T5),mean(T6),mean(T7),mean(T8),mean(T9),mean(T10)]
stdT = [std(T1),std(T2),std(T_alpha),std(T_beta),std(T3),std(T4),std(T5),std(T6),std(T7),std(T8),std(T9),std(T10)]
end
