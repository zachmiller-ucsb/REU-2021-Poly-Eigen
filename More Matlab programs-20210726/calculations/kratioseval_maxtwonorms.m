function maxB = kratioseval_maxtwonorms(bks)
% Return the maximum two norm of the matrix coefficients bks.

    bklength = size(bks, 3);
    twonormbounds = zeros(1, bklength);
    for k = 1:bklength
        twonormbounds(k) = stdtwonormmax(bks(:,:,k));
    end
    [twonormbounds, twonormorder] = sort(twonormbounds, 'descend');
    maxB = stdtwonorm(bks(:,:,twonormorder(1)));
    for k = 2:bklength
        if maxB < twonormbounds(twonormorder(k))
            maxB = max(maxB, stdtwonorm(bks(:,:,twonormorder(k))));
        else
            break;
        end
    end
end
