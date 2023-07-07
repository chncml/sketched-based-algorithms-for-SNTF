function [G, Obj, t, fit] = uniformrand_Parallel_Multi_SNTF(Y, G0, alpha, maxIter, maxTime, tol, tau, r, computeobj)

istol = 0;
if tol ~= 0
   istol = 1;
end

[nrow, ncol] = size(G0);snrow = nrow^2;
Ymat = reshape(Y, [nrow, snrow]);
prob = ones(snrow, 1)/snrow;
samplerow = ceil(tau*snrow);
idx = randsample(snrow, samplerow, true, prob);
Ymat1 = Ymat(:, idx)/sqrt(snrow/samplerow);
omega = rand(samplerow, ncol + r)/sqrt(ncol + r);
Ymat1 = Ymat1 * omega;

G = G0;
for iter = 1:maxIter    

    tstart = tic;
      
    % compute columnwise Kronecker product
    khatriraoG = kr(G, G);
    khatriraoG1 = omega' * khatriraoG(idx, :)/sqrt(snrow/samplerow);
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
