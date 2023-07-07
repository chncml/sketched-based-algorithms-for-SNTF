function [G, Obj, t, fit] = Parallel_Multi_SNTF(Y, G0, alpha, maxIter, maxTime, tol, computeobj)

istol = 0;
if tol ~= 0
   istol = 1;
end

[nrow, ncol] = size(G0);
Ymat = reshape(Y,[nrow,nrow^2]); 

G = G0;
khatriraoG = zeros(nrow^2, ncol);
for iter = 1:maxIter    

    tstart = tic;
      
    GtG = G'*G; 
    % compute columnwise Kronecker product
    khatriraoG = kr(G, G);
    numerator = Ymat * khatriraoG;
    % numerator = Ymat * khatrirao(G,G);
    denominator = G*(GtG.*GtG); 
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
