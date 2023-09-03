function [A_new, b_vec_new] = remove_null_row(A,b_vec)
    for i = 1:size(A,1)
        if norm(A(i,:))<1e-3 && b_vec(i)<-1e-3
            error('The %th row is problamatic',i)
        end
    end
    row_A =vecnorm(A,2,2)<1e-3.*ones(size(A,1),1);
    row_b =abs(b_vec) > 1e-4;
    row_to_be_removed = row_A;
   
    A(row_to_be_removed,:) = [];
    b_vec(row_to_be_removed,:) = [];
    
    A_new = A;
    b_vec_new = b_vec;
end
