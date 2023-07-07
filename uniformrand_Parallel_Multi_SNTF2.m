function [G, Obj, t, fit] = uniformrand_Parallel_Multi_SNTF2(Y, G0, alpha, maxIter, maxTime, tol, tau, r, computeobj)

istol = 0;
if tol ~= 0
   istol = 1;
end

[nrow, ncol] = size(G0);
prob = ones(nrow, 1)/nrow;
samplerow = ceil(tau * nrow);
ind = randsample(nrow, samplerow, true, prob);
Y1 = Y(:, ind, ind)/(nrow/samplerow);
Ymat = reshape(Y1, nrow, []);
omega = rand(samplerow, ncol + r)/sqrt(ncol + r);
Ymat1 = Ymat * kr(omega, omega);

G = G0;
for iter = 1:maxIter    

    tstart = tic;
      
    % compute columnwise Kronecker product
    G1 = G(ind, :)/sqrt(nrow/samplerow);
    GtG = omega'  * G1;
    khatriraoG1 = GtG.*GtG;
    numerator = Ymat1 * khatriraoG1;
    % numerator = Ymat * khatrirao(G,G);
    denominator = G * (khatriraoG1' * khatriraoG1); 
    G = G.* ((numerator./denominator).^(alpha)); 
    G = max(G,1e-16);
    
    t(iter) = toc(tstart);
        
    % record the objective function value of SNCP
    if  computeobj % && mod(iter,10)==1
        
        Ztensor = ktensor({G,G,G}); 
        Z = double(tensor(Ztensor));
        Obj(iter) = norm(Y(:)-Z(:))^2;      

    end

    tstart = tic;
    
    if istol && norm(G-G0,'fro')<tol*norm(G0,'fro')
       break;
    end
        
    if sum(t)>maxTime
       break;
    end
    
    G0 = G;
    
    t(iter) = t(iter) + toc(tstart);
    
end
t = cumsum(t);

% compute the fit degree
Ztensor = ktensor({G,G,G}); 
Z = double(tensor(Ztensor));
fit = 1 - norm( Y(:)-Z(:) ) / norm( Y(:) );

if ~computeobj
   Obj = norm(Y(:)-Z(:))^2;
end
