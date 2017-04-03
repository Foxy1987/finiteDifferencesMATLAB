function [ r ] = tridiag( a, b, c, r)
%%Solve equation Qu^{n+1} =  by writing the augmented matrix [Q|r] and use
%%Gaussian Elimination. Write tridiagonal matrix in the form T[a, b, c]
%%with the terms a_2, ..., a_m on the subdiagonal, the terms b_1, ..., b_m
%%on the diagonal and the terms c_1, ..., c_{m-1} on the superdiagonal,
%%then use the Thomas algorithm to solve the system
%%
    m=length(r);
    c(1) = c(1)/b(1);
    r(1) = r(1)/b(1);
    
    for j = 2:m-1
        temp = b(j) - a(j)*c(j-1);
        r(j) = r(j) - a(j)*r(j-1);
        c(j) = c(j)/temp;
        r(j) = r(j)/temp;
    end
    
    b(m) = b(m) - a(m)*c(m-1);
    r(m) = r(m) - a(m)*r(m-1);
    r(m) = r(m)/b(m);
    
    for j = m-1:-1:1
       r(j) = r(j) - c(j)*r(j+1); 
    end
    
end

