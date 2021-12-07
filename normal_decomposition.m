function [Ynext,Xcurr,P,cnt] = normal_decomposition(X,Ycurr,L,N,Niter)

    Xcurr = [];
    Mcurr = zeros(2,L);
    Sigma_curr_vector = zeros(2,2,L);
    P = zeros(1,L);

    for i = 1:L
        Xcurr_i = X(1:2,Ycurr==i);
        Xcurr = [Xcurr, Xcurr_i];
        Mcurr(:,i) = mean(Xcurr(i),2);
        Sigma_curr_vector(:,:,i) = cov(Xcurr_i')';
        P(i) = length(Xcurr_i)/length(X);
    end
    
    
    condition = 0;
    Ynext = zeros(1,L*N);
    cnt = 0;
    while (condition == 0 && Niter == 0) || (condition == 0 && cnt < Niter && Niter > 0)
        cnt = cnt+1;
        for j = 1:2*N
            J = zeros(1,L);
            for i = 1:L
                delta = X(1:2,j) - Mcurr(:,i);
                Sigma_i = Sigma_curr_vector(:,:,i);
                J(i) = -2*log(P(i)) + log(det(Sigma_i)) + delta'*Sigma_i^(-1)*delta;
            end
            [~, Ynext(j)] = min(J);
        end

        Xcurr = [];
        Mcurr = zeros(2,L);    
        Sigma_curr_vector = zeros(2,2,L);

        for i = 1:L
            Xcurr_i = X(1:2,Ynext==i);
            Mcurr(:,i) = mean(Xcurr_i,2);
            Xcurr = [Xcurr, Xcurr_i];
            Sigma_curr_vector(:,:,i) = cov(Xcurr_i')';
            P(i) = length(Xcurr_i)/length(X);
        end

        condition = isequal(Ycurr,Ynext);
        Ycurr = Ynext;
        

    end

    if (L == 2)
        N1 = 1;
        N2 = 1000;
        Ynext(Ycurr==Ycurr(N1)) = 1;
        Ynext(Ycurr==Ycurr(N2)) = 2;

        P = sum(X(3,:)==Ynext)/2/N;
    else
        P = 1;
    end



end