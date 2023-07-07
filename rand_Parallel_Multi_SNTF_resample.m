function [G, Obj, t, fit] = rand_Parallel_Multi_SNTF_resample(Y, G0, alpha, maxIter, maxTime, tol, r, computeobj)
istol = 0;
if tol ~= 0
   istol = 1;
end
[nrow, ncol] = size(G0);snrow = nrow^2;
Ymat1 = reshape(Y,[nrow,snrow]);
G = G0;
for iter = 1:maxIter    
    tstart = tic;
    omega = rand(snrow, ncol + r)/sqrt(ncol + r);
    Ymat = Ymat1 * omega;
    G1 = kr(G, G);GtG = omega' * G1; 
    % compute columnwise Kronecker product
    numerator = Ymat * GtG;
    % numerator = Ymat * khatrirao(G,G);
    denominator = G*(GtG'*GtG); 
    G = G.* (numerator./denominator).^(alpha); 
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
