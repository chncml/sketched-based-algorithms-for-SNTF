function [G,f,t,fit] = beta_Paralle_Multi_SNTF(Y,G0,beta,maxIter,maxTime,tol)

tstart = tic;

istol = 0;
if tol ~= 0
   istol = 1;
end

[nrow, ncol] = size(G0);
Ymat = reshape(Y,[nrow,nrow^2]); 
normY = norm(Y(:))^2;

% compute the objective value of the first iteration
G = G0;  GtG = G'*G; 
khatriraoG = zeros(nrow^2,ncol);
for i = 1:ncol
    temp = G(:,i)*G(:,i)';
    khatriraoG(:,i) = temp(:);
end
numerator = Ymat * khatriraoG;
% numerator = Ymat * khatrirao(G,G);
denominator = G*(GtG.*GtG); 
vecG = G(:)';
Obj0 = normY - 2*(vecG * numerator(:)) + vecG * denominator(:);

d = 1-beta; 
for iter = 1:maxIter    
    
    R = numerator./denominator;
    Gnew = G.* (d+beta*R); 
    Gnew = max(Gnew,1e-16);
    
    % compute the objective value
    GtG = Gnew'*Gnew;
    for i = 1:ncol
        temp = Gnew(:,i)*Gnew(:,i)';
        khatriraoG(:,i) = temp(:);
    end
    numerator = Ymat * khatriraoG;
    denominator = Gnew*(GtG.*GtG); 
    vecG = Gnew(:)';
    Obj = normY - 2*(vecG * numerator(:)) + vecG * denominator(:);
    
    if Obj >= Obj0
       
       Gnew = G.* (R.^(1/5));
       Gnew = max(Gnew,1e-16);
       
       % compute the objective value
       GtG = Gnew'*Gnew; 
       for i = 1:ncol
           temp = Gnew(:,i)*Gnew(:,i)';
           khatriraoG(:,i) = temp(:);
       end
       numerator = Ymat * khatriraoG;
       denominator = Gnew*(GtG.*GtG); 
       vecG = Gnew(:)';
       Obj = normY - 2*(vecG * numerator(:)) + vecG * denominator(:);
       
    end   
   
    f(iter) = Obj;
    t(iter) = toc(tstart);
    
    if istol && norm(Gnew-G,'fro')<tol*norm(G,'fro')
       break;
    end
   
    if t(end)>maxTime
       break;
    end
    
    G = Gnew;
    Obj0 = Obj;   
    
end
G = Gnew; 

% compute the fit degree
Ztensor = ktensor({G,G,G}); 
Z = double(tensor(Ztensor));
fit = 1 - norm( Y(:)-Z(:) ) / norm( Y(:) );


