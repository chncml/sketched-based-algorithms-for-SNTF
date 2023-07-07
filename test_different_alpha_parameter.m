function test_different_alpha_parameter
clc; clear; close all
format short e
I = 200;J = 20;R = J;SNR = 10; 

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
    tic;[G1, f1, t1, fit1] = rand_Parallel_Multi_SNTF(Y, G0, (1/6), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G2, f2, t2, fit2] = rand_Parallel_Multi_SNTF(Y, G0, (1/5), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G3, f3, t3, fit3] = rand_Parallel_Multi_SNTF(Y, G0, (1/4), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G4, f4, t4, fit4] = rand_Parallel_Multi_SNTF(Y, G0, (1/3), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G5, f5, t5, fit5] = rand_Parallel_Multi_SNTF(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G6, f6, t6, fit6] = rand_Parallel_Multi_SNTF(Y, G0, (3/5), maxIter, maxTime, tol, 20, computeobj);toc;
    
    T11(MCT,1) = t1(end);
    T12(MCT,1) = t2(end);
    T13(MCT,1) = t3(end);
    T14(MCT,1) = t4(end);
    T15(MCT,1) = t5(end);
    T16(MCT,1) = t6(end);

    Fit11(MCT,1) = fit1(end);
    Fit12(MCT,1) = fit2(end);
    Fit13(MCT,1) = fit3(end);
    Fit14(MCT,1) = fit4(end);
    Fit15(MCT,1) = fit5(end);
    Fit16(MCT,1) = fit6(end);

    MCT = MCT+1;

end
 
% Monte Carlo tests
MCT = 1; NumofMCT = 10; 
while MCT <= NumofMCT
    MCT
    
        % initialization
     G0 = rand(I, R) + 1e-5;
%     G0 = ones(I, R);
  
    tic;[G1, f1, t1, fit1] = randkr_Parallel_Multi_SNTF1(Y, G0, (1/6), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G2, f2, t2, fit2] = randkr_Parallel_Multi_SNTF1(Y, G0, (1/5), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G3, f3, t3, fit3] = randkr_Parallel_Multi_SNTF1(Y, G0, (1/4), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G4, f4, t4, fit4] = randkr_Parallel_Multi_SNTF1(Y, G0, (1/3), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G5, f5, t5, fit5] = randkr_Parallel_Multi_SNTF1(Y, G0, (1/2), maxIter, maxTime, tol, 20, computeobj);toc;
    tic;[G6, f6, t6, fit6] = randkr_Parallel_Multi_SNTF1(Y, G0, (3/5), maxIter, maxTime, tol, 20, computeobj);toc;
    
    T21(MCT,1) = t1(end);
    T22(MCT,1) = t2(end);
    T23(MCT,1) = t3(end);
    T24(MCT,1) = t4(end);
    T25(MCT,1) = t5(end);
    T26(MCT,1) = t6(end);

    Fit21(MCT,1) = fit1(end);
    Fit22(MCT,1) = fit2(end);
    Fit23(MCT,1) = fit3(end);
    Fit24(MCT,1) = fit4(end);
    Fit25(MCT,1) = fit5(end);
    Fit26(MCT,1) = fit6(end);

    MCT = MCT+1;

end
        
disp('The results are shown as follows:')


meanFIT1 = [mean(Fit11),mean(Fit12),mean(Fit13),mean(Fit14),mean(Fit15),mean(Fit16)]
stdFIT1 = [std(Fit11),std(Fit12),std(Fit13),std(Fit14),std(Fit15),std(Fit16)]

meanT1 = [mean(T11),mean(T12),mean(T13),mean(T14),mean(T15),mean(T16)]
stdT1 = [std(T11),std(T12),std(T13),std(T14),std(T15),std(T16)]

meanFIT2 = [mean(Fit21),mean(Fit22),mean(Fit23),mean(Fit24),mean(Fit25),mean(Fit26)]
stdFIT2 = [std(Fit21),std(Fit22),std(Fit23),std(Fit24),std(Fit25),std(Fit26)]

meanT2 = [mean(T21),mean(T22),mean(T23),mean(T24),mean(T25),mean(T26)]
stdT2 = [std(T21),std(T22),std(T23),std(T24),std(T25),std(T26)]
end
