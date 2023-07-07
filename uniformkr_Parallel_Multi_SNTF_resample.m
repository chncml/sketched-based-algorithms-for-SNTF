function [G, Obj, t, fit] = uniformkr_Parallel_Multi_SNTF_resample(Y, G0, alpha, maxIter, maxTime, tol, tau, computeobj)

istol = 0;
if tol ~= 0
   istol = 1;
end

[nrow, ~] = size(G0);
% prob = sum(Ymat.^2, 2)/sum(Ymat(:).^2);
prob = ones(nrow, 1)/nrow;
samplerow = ceil(tau * nrow);


G = G0;
for iter = 1:maxIter    

    tstart = tic;
    ind = randsample(nrow, samplerow, true, prob);
    Y1 = Y(:, ind, ind)/(nrow/samplerow);
    Ymat = reshape(Y1, nrow, []);
    G1 = G(ind, :)/sqrt(nrow/samplerow);
    GtG = kr(G1, G1); 
    numerator = Ymat * GtG;denominator = G*(GtG'*GtG); 
    G = G.* ((numerator./denominator).^(alpha)); 
    G = max(G,1e-16);t(iter) = toc(tstart);
        
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
