function [ r ] = tridiag2( a, b, c, r, n)
%%Solve equation Qu^{n+1} =  by writing the augmented matrix [Q|r] and use
%%Gaussian Elimination. Write tridiagonal matrix in the form T[a, b, c]
%%with the terms a_2, ..., a_m on the subdiagonal, the terms b_1, ..., b_m
%%on the diagonal and the terms c_1, ..., c_{m-1} on the superdiagonal,
%%then use the Thomas algorithm to solve the system

%%n denotes the number of RHS equations

%%
    m=length(r);
    c(1) = c(1)/b(1);
    
    % each column is a RHS vector, loop through
    for i =1:n
        r(1,i) = r(1,i)/b(1);
    end
    
    % only need to index the RHS columns because RHS vectors share the same
    % coefficient matrix
    for i =1:n
        for j = 2:m-1
            temp = b(j) - a(j)*c(j-1);
            r(j, i) = r(j, i) - a(j)*r(j-1, i);
            c(j) = c(j)/temp;
            r(j, i) = r(j, i)/temp;
        end   
    end
    
    b(m) = b(m) - a(m)*c(m-1);
    for i=1:n
        r(m,i) = r(m,i) - a(m)*r(m-1,i);
        r(m,i) = r(m,i)/b(m);
    end
    
    for i=1:n
        for j = m-1:-1:1
            r(j,i) = r(j,i) - c(j)*r(j+1,i); 
        end
    end
    
end
