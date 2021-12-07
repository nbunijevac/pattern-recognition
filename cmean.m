function [Ynext,Xcurr,P,cnt] = cmean(X,Ycurr,L,N,Niter)

    Xcurr = [];
    Mcurr = zeros(2,L);

    for i = 1:L
        Xcurr_i = X(1:2,Ycurr==i);
        Xcurr = [Xcurr, Xcurr_i];
        Mcurr(:,i) = mean(Xcurr(i),2);
    end

    condition = 0;
    Ynext = zeros(1,L*N);
    cnt = 0;
    while (condition == 0 && Niter == 0) || (condition == 0 && cnt < Niter && Niter > 0)
        cnt = cnt+1;
        for i = 1:4*N
            norm_vector = zeros(1,L);
            for j = 1:L
                norm_vector(j) = norm(X(1:2,i) - Mcurr(:,j));
            end
            [~, Ynext(i)] = min(norm_vector);
        end

        Xcurr = [];
        Mcurr = zeros(2,L);
        for i = 1:L
            Xcurr_i = X(1:2,Ynext==i);
            Mcurr(:,i) = mean(Xcurr_i,2);
            Xcurr = [Xcurr, Xcurr_i];
        end

        condition = isequal(Ycurr,Ynext);
        Ycurr = Ynext;

    end

    if (L == 4)
        N1 = 144;
        N2 = 525;
        N3 = 1022;
        N4 = 1652;
        Ynext(Ycurr==Ycurr(N1)) = 1;
        Ynext(Ycurr==Ycurr(N2)) = 2;
        Ynext(Ycurr==Ycurr(N3)) = 3;
        Ynext(Ycurr==Ycurr(N4)) = 4;

        P = sum(X(3,:)==Ynext)/4/N;
    else
        P = 1;
    end



end