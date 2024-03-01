function M = GetMarkovMatrix(C, A, B, pp)
%GETMARKOVMATRIX returns a matrix of structure pp such that
% M = [CA^pp(1,1)B, CA^pp(1,2)B ... CA^pp(1,end)B;
%      :                :                   :    ;
%      CA^pp(end,1)B, ...       ... CA^pp(end,end)B];
%
% Further exceptions: if any value of pp < 0, the corresponding block will
% be zero instead. use pp = 0 for CB

    m = size(C,1);
    n = size(B,2);

    M = zeros(size(pp,1)*m, size(pp,2)*n);

    for i = 1:size(pp,1)
        for j = 1:size(pp,2)
            if(pp(i,j) >= 0)
                M((i-1)*m+1:i*m, (j-1)*n+1:j*n) = C*A^pp(i,j)*B;
            end
        end
    end

end

