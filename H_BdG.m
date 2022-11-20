function op = H_BdG(N_row,N_col,Jx,Jy,Jz,K,t,Del_ep,mate,z1,z2)

% we start constructing the Hamiltonian with the Majorana representations,
% and finally transform to that of the fermions.

% we consider the type of bonds one by one, and we first construct the
% A matrix in related to the Majorana operators: H = I*A \gammga \gamma'

% --- x-bonds --- %

Ax = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Ax( findex(x,y,N_col), findex(fhop(x+1,N_row),y,N_col) ) = -Jx*(-1)^sum( mate(x,y+1:N_col) )...
            *z2^(x == N_row);
    end
%  right boundary no e parity term
    Ax(findex(x,N_col,N_col), findex(fhop(x+1,N_row),N_col,N_col)) = -Jx*z2^(x == N_row);
end

% --- y-bonds --- %

Ay = zeros(N_row*N_col);

for x = 1:N_row-1
    for y = 1:N_col-1
        Ay(findex(x,y,N_col), findex(x+1,y+1,N_col)) = Jy*(-1)*(-1)^(sum(mate(x,y+1:N_col)));
    end
    Ay( findex(x,N_col,N_col), findex(x+1,1,N_col)) = Jy*(-1)*(-1)^( sum(mate(x+1:N_row,:), [1 2] ) )...
        *z1;
end

for y = 1:N_col-1
    Ay(findex(N_row,y,N_col), findex(1,y+1,N_col)) = Jy*(-1)*(-1)^(sum(mate(N_row,y+1:N_col)))*z2;
end

Ay(findex(N_row,N_col,N_col), findex(1,1,N_col)) = Jy*(-1)*z1*z2;

% --- z-bonds --- %

Az = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Az( findex(x,y,N_col), findex(x,y+1,N_col) ) = -Jz;
    end
    Az( findex(x,N_col,N_col), findex(x,1,N_col) ) = -Jz*(-1)^( sum(mate(x:N_row,:), [1 2] ) )*z1;
end


% --- the 3rd neighbour spin couplings --- %
% remember that we need to time a "-t/2" finally

% --- horizontal direction --- %

At1 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-2
%         At1( findex(x,y,N_col), findex(x,y+2,N_col) ) = 1+1;
        At1( findex(x,y,N_col), findex(x,y+2,N_col) ) = 1 + (-1)^mate(x,y+1);
    end
%     At1( findex(x,N_col-1,N_col), findex(x,1,N_col) ) = ( 1+1 )*...
%         (-1)^( sum(mate(x:N_row,:), [1 2]) )*z1;
%     At1( findex(x,N_col,N_col), findex(x,2,N_col) ) = ( 1+1 )*...
%         (-1)^( sum(mate(x:N_row,:), [1 2]) )*z1;
    At1( findex(x,N_col-1,N_col), findex(x,1,N_col) ) = ( 1+(-1)^mate(x,N_col) )*...
        (-1)^( sum(mate(x:N_row,:), [1 2]) )*z1;
    At1( findex(x,N_col,N_col), findex(x,2,N_col) ) = ( 1+(-1)^mate(x,1) )*...
        (-1)^( sum(mate(x:N_row,:), [1 2]) )*z1;
end

% --- vertical direction --- %

At2 = zeros(N_row*N_col);

for x = 1:N_row-2
    for y = 1:N_col-1
%         At2( findex(x,y,N_col), findex(fhop(x+2,N_row),y,N_col) ) = ( 1+1 )*...
%             (-1)^( sum(mate(x,y+1:N_col)) + sum(mate(fhop(x+1,N_row),y+1:N_col)) );
        At2( findex(x,y,N_col), findex(fhop(x+2,N_row),y,N_col) ) = ( 1+(-1)^mate(fhop(x+1,N_row),y) )*...
            (-1)^( sum(mate(x,y+1:N_col)) + sum(mate(fhop(x+1,N_row),y+1:N_col)) );
    end
%     At2( findex(x,N_col,N_col), findex(fhop(x+2,N_row),N_col,N_col) ) = ( 1+1 );
    At2( findex(x,N_col,N_col), findex(fhop(x+2,N_row),N_col,N_col) ) = ( 1+(-1)^mate(fhop(x+1,N_row),N_col) );
end

for x = N_row-1:N_row
    for y = 1:N_col-1
%         At2( findex(x,y,N_col), findex(fhop(x+2,N_row),y,N_col) ) = ( 1+1 )*...
%             (-1)^( sum(mate(x,y+1:N_col)) + sum(mate(fhop(x+1,N_row),y+1:N_col)) )*z2;
        At2( findex(x,y,N_col), findex(fhop(x+2,N_row),y,N_col) ) = ( 1+(-1)^mate(fhop(x+1,N_row),y) )*...
            (-1)^( sum(mate(x,y+1:N_col)) + sum(mate(fhop(x+1,N_row),y+1:N_col)) )*z2;
    end
%     At2( findex(x,N_col,N_col),findex(fhop(x+2,N_row),N_col,N_col) ) = ( 1+1 )*z2;
    At2( findex(x,N_col,N_col),findex(fhop(x+2,N_row),N_col,N_col) ) = ( 1+(-1)^mate(fhop(x+1,N_row),N_col) )*z2;
end

% --- the coupling within a plaquette --- %

At3 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col
%         At3( findex(x,y,N_col), findex(x,y,N_col) ) = 1+1;
        At3( findex(x,y,N_col), findex(x,y,N_col) ) = 1+(-1)^mate(x,y);
    end
end

% --- the Del_ep is a local potential term Del_ep c^\dagger c, so the
% Del_ep is essentially the (-1)*\mu, which \mu being the chemical
% potential --- %

A = Ax + Ay + Az + Del_ep/2*eye(N_row*N_col) + (-t/2)*( At1 + At2 + At3 );


% --- the D1 and D2 defined below is for the case with some local chemical
% potential around the added e particles, only use this for the C = -2
% case, for the C = 1 or pure Kitaev model, do not use it and just use the
% one above --- %

% D1 = zeros(N_row*N_col);
% D2 = zeros(N_row*N_col);
% 
% [ x1, y1 ] = find( mate == 1 );
% [ x2, y2 ] = find( mate == -1 );
% 
% if length(x1) == 1
%     D1(findex(x1,y1,N_col), findex(x1,y1,N_col)) = Del_ep/2;
%     D1(findex(fhop(x1+1,N_row),y1,N_col), findex(fhop(x1+1,N_row),y1,N_col)) = Del_ep/2;
%     D1(findex(x1,fhop(y1-1,N_col),N_col), findex(x1,fhop(y1-1,N_col),N_col)) = Del_ep/2;
%     D1(findex(fhop(x1+1,N_row),fhop(y1-1,N_col),N_col), findex(fhop(x1+1,N_row),fhop(y1-1,N_col),N_col)) = Del_ep/2;
% %     
%     D2(findex(x2,y2,N_col), findex(x2,y2,N_col)) = -Del_ep/2;
%     D2(findex(fhop(x2+1,N_row),y2,N_col), findex(fhop(x2+1,N_row),y2,N_col)) = -Del_ep/2;
%     D2(findex(x2,fhop(y2-1,N_col),N_col), findex(x2,fhop(y2-1,N_col),N_col)) = -Del_ep/2;
%     D2(findex(fhop(x2+1,N_row),fhop(y2-1,N_col),N_col), findex(fhop(x2+1,N_row),fhop(y2-1,N_col),N_col)) = -Del_ep/2;
%     % D1( findex(fix(N_row/2)+1,1,N_col), findex(fix(N_row/2)+1,1,N_col) ) = Del_ep/2*mate(fix(N_row/2)+1, 1);
%     % D1( findex(fix(N_row/2)+2,1,N_col), findex(fix(N_row/2)+2,1,N_col) ) = Del_ep/2*mate(fix(N_row/2)+1, 1);
%     % D1( findex(fix(N_row/2)+1,N_col,N_col), findex(fix(N_row/2)+1,N_col,N_col) ) = Del_ep/2*mate(fix(N_row/2)+1, 1);
%     % D1( findex(fix(N_row/2)+2,N_col,N_col), findex(fix(N_row/2)+2,N_col,N_col) ) = Del_ep/2*mate(fix(N_row/2)+1, 1);
%     % 
%     % for j = 2:N_col
%     %     D2( findex(fix(N_row/2)+1,j,N_col), findex(fix(N_row/2)+1,j,N_col) ) = -Del_ep/2*mate(fix(N_row/2)+1, j);
%     %     D2( findex(fix(N_row/2)+2,j,N_col), findex(fix(N_row/2)+2,j,N_col) ) = -Del_ep/2*mate(fix(N_row/2)+1, j);
%     %     D2( findex(fix(N_row/2)+1,j-1,N_col), findex(fix(N_row/2)+1,j-1,N_col) ) = -Del_ep/2*mate(fix(N_row/2)+1, j);
%     %     D2( findex(fix(N_row/2)+2,j-1,N_col), findex(fix(N_row/2)+2,j-1,N_col) ) = -Del_ep/2*mate(fix(N_row/2)+1, j);
%     % end
%     D = D1 + D2;
% else
%     D = D1 + D2;
% end

% A = Ax + Ay + Az + D + (-t/2)*( At1 + At2 + At3 );

% ---end--- %


% note that the A here is not the A matrix in Kitaev's paper!
% see the end!


% the couplings between the i gamma & gammas:
% remember there is a "-K" factor before, which we will multiply at the end

% horizontal couplings:

Ag1 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Ag1(findex(x,y,N_col),findex(x,y+1,N_col)) = (-1)^mate(x,y+1);
    end
    Ag1(findex(x,N_col,N_col),findex(x,1,N_col)) =...
        (-1)^( sum(mate(x:N_row,:), [1 2] )+mate(x,1) )*z1;
end

% verticle couplings

Ag2 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Ag2(findex(x,y,N_col),findex(fhop(x+1,N_row),y,N_col)) =...
            (-1)*(-1)^( sum( mate(x,y+1:N_col) ) )*z2^(x == N_row);
    end
    Ag2(findex(x,N_col,N_col),findex(fhop(x+1,N_row),N_col,N_col)) = (-1)*z2^(x == N_row);
end

% off-diagonal sites

Ag3 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Ag3( findex(x,y,N_col),findex(fhop(x+1,N_row),fhop(y-1,N_col),N_col) ) =...
            (-1)^( sum(mate(x,y+1:N_col)) + (y==1)*sum(mate(fhop(x+1,N_row):N_row,:), [1 2] ) )...
            *z1^(y==1)*z2^(x==N_row);
    end
    Ag3(findex(x,N_col,N_col),findex(fhop(x+1,N_row),N_col-1,N_col)) = z2^(x==N_row);
end

Ag = -K*( (Ag1 + Ag2 + Ag3) - (Ag1 + Ag2 + Ag3).' );
% so in the Hamiltonina, this part contribute as: i/2 gammga*Ag*gamma.

% the gamma'&gamma' coupling:

% horizontal couplings:

Agp1 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Agp1(findex(x,y,N_col),findex(x,y+1,N_col)) = -1;
    end
    Agp1(findex(x,N_col,N_col),findex(x,1,N_col)) = -1*(-1)^(sum(mate(x:N_row,:), [1 2] ))*z1;
end

% verticle couplings:

Agp2 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col
        Agp2(findex(x,y,N_col),findex(fhop(x+1,N_row),y,N_col)) =...
            (-1)^( sum(mate(x,y:N_col)) )*z2^( x == N_row );
    end
end

% off-diagonal coupling:

Agp3 = zeros(N_row*N_col);

for x = 1:N_row
    for y = 1:N_col-1
        Agp3(findex(x,y,N_col),findex(fhop(x-1,N_row),y+1,N_col)) =...
            (-1)^( sum(mate(fhop(x-1,N_row),y+1:N_col)) )*z2^(x==1);
    end
    Agp3(findex(x,N_col,N_col),findex(fhop(x-1,N_row),1,N_col)) =...
        (-1)^( sum(mate(fhop(x-1,N_row):N_row, :), [1 2] ) )*z1*z2^(x==1);
end

Agp = -K*( (Agp1 + Agp2 + Agp3) - (Agp1 + Agp2 + Agp3).' );

% now we have: H = i/2 (\gamma, \gamma' ) ( Ag, A; -A^T, Agp ) ( \gamma; \gamma' )
% the BdG Hamiltonian

op = 1i*[ eye(N_row*N_col), 1i*eye(N_row*N_col); eye(N_row*N_col), -1i*eye(N_row*N_col) ]*...
    [ Ag, A; -A.', Agp ]*...
    [ eye(N_row*N_col), eye(N_row*N_col); -1i*eye(N_row*N_col), 1i*eye(N_row*N_col) ];


end