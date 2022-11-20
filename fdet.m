function op = fdet(A)


[ L, U ] = lu(A);

s =  det(L);


% --- calcuating the product by taking the exponential of the sum
% of logarithmic of each elements, and taking sum will be free of the
% blowing up issue --- %


nada = log( diag(U) );

op = s*exp( sum(nada) );


% ---end--- %



% --- product of the descended and ascended list to avoid blowing up/vanishing
% of the product --- %


% nada = sort(diag(U), "ascend") .* sort(diag(U), "descend");
% 
% op = s*sqrt( prod(nada) );


% --- end --- %



end