function test_different_RS_probability_model
clc; clear; close all
format short e
tau = 0.5;

I = 400;J = 20;R = J;SNR = 0; 

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
    tic;[G1, f1, t1, fit1] = uniform_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, tau^2, computeobj);toc;
    tic;[G2, f2, t2, fit2] = uniformkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, tau, computeobj);toc;
    tic;[G3, f3, t3, fit3] = uniform_Parallel_Multi_SNTF_optimal(Y, G0, (1/2), maxIter, maxTime, tol, tau^2, computeobj);toc;
    tic;[G4, f4, t4, fit4] = uniformkr_Parallel_Multi_SNTF_optimal(Y, G0, (1/2), maxIter, maxTime, tol, tau, computeobj);toc;
    
    T11(MCT,1) = t1(end);
    T12(MCT,1) = t2(end);
    T13(MCT,1) = t3(end);
    T14(MCT,1) = t4(end);

    Fit11(MCT,1) = fit1(end);
    Fit12(MCT,1) = fit2(end);
    Fit13(MCT,1) = fit3(end);
    Fit14(MCT,1) = fit4(end);

    MCT = MCT+1;

end

meanFIT1 = [mean(Fit11), mean(Fit12), mean(Fit13), mean(Fit14)]
stdFIT1 = [std(Fit11), std(Fit12), std(Fit13), std(Fit14)]

meanT1 = [mean(T11), mean(T12), mean(T13), mean(T14)]
stdT1 = [std(T11), std(T12), std(T13), std(T14)]

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
    tic;[G1, f1, t1, fit1] = uniform_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, tau^2, computeobj);toc;
    tic;[G2, f2, t2, fit2] = uniformkr_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, tau, computeobj);toc;
    tic;[G3, f3, t3, fit3] = uniform_Parallel_Multi_SNTF_optimal(Y, G0, (1/2), maxIter, maxTime, tol, tau^2, computeobj);toc;
    tic;[G4, f4, t4, fit4] = uniformkr_Parallel_Multi_SNTF_optimal(Y, G0, (1/2), maxIter, maxTime, tol, tau, computeobj);toc;
    
    T11(MCT,1) = t1(end);
    T12(MCT,1) = t2(end);
    T13(MCT,1) = t3(end);
    T14(MCT,1) = t4(end);

    Fit11(MCT,1) = fit1(end);
    Fit12(MCT,1) = fit2(end);
    Fit13(MCT,1) = fit3(end);
    Fit14(MCT,1) = fit4(end);

    MCT = MCT+1;

end

meanFIT1 = [mean(Fit11), mean(Fit12), mean(Fit13), mean(Fit14)]
stdFIT1 = [std(Fit11), std(Fit12), std(Fit13), std(Fit14)]

meanT1 = [mean(T11), mean(T12), mean(T13), mean(T14)]
stdT1 = [std(T11), std(T12), std(T13), std(T14)]
end
