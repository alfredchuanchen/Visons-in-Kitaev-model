load('parameters.mat');

sp = 1;

% lambda = 2.85 - 1.85*tanh(abs(Del_ep)/1.1); % for AFM case with small K
% lambda = 2.21 - 1.21*tanh(abs(Del_ep)/1.3);

% lambda = 1.61 - 0.61*tanh(abs(Del_ep)/2); % for FM case
lambda = 1.61;

% this is a parameter used to increase the determinant and the pfaffian 
% to to avoid the numerical 0.


% ---when moving one particle, the row and column coordinate along the path---%

x_i = path( sp, 1 ); % row
y_i = path( sp, 2 ); % column
x_f = path( mod(sp,L)+1, 1 );
y_f = path( mod(sp,L)+1, 2 );


% ---braiding of the two particle in a gauge invariant way--- %

% if sp < L
%     nada = abs( path(sp+1,:) - path(sp,:) );
%     loc = find(nada == 1);
%     if mod(loc,2) == 1
%         x_i = path( sp, loc ); % row
%         y_i = path( sp, loc+1 ); % column
%         x_f = path( sp+1, loc );
%         y_f = path( sp+1, loc+1 );
%     else
%         x_i = path( sp, loc-1);
%         y_i = path( sp, loc );
%         x_f = path( sp+1, loc-1 );
%         y_f = path( sp+1, loc );
%     end
% else
% % ---the final step for the exchange of the 2 e particles--- %
% %     x_i = path( sp, 1 ); % row
% %     y_i = path( sp, 2 ); % column
% %     x_f = path( 1, 3 );
% %     y_f = path( 1, 4 );
% % ---end--- %
% % 
% % ---the final step for the double exchange--- %
%     x_i = path( sp, 3 ); % row
%     y_i = path( sp, 4 ); % column
%     x_f = path( 1, 3 );
%     y_f = path( 1, 4 );
% % ---end--- %
% end

% this is the "bra" state (the state after acting the Y operaotr)

load( strcat( "spectra-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string( mod(sp,L)+1 ),...
    ".mat" ) );
% ---notice here the directory of the eigenvector data: in the cluster, it
% is stored in the same directory as this script

mat1 = poseigenvec;
u1 = mat1( 1:N_row*N_col, 1:N_row*N_col );
v1 = mat1( N_row*N_col+1 : 2*N_row*N_col, 1:N_row*N_col );
clear poseigenvec;

% this is the ket state (the state before acting the Y operator)

load( strcat( "spectra-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string( sp ),...
    ".mat" ) );

mat2 = poseigenvec;
u2 = mat2( 1:N_row*N_col, 1:N_row*N_col );
v2 = mat2( N_row*N_col+1 : 2*N_row*N_col, 1:N_row*N_col );
clear poseigenvec;

I_p = eye(N_row*N_col);
for y = 1 : min( [y_i, y_f] ) - 1
    I_p( findex( max([ x_i, x_f ]), y, N_col ), findex( max([ x_i, x_f ]), y, N_col ) ) = -1;
end
newu2 = I_p*u2*I_p;
newv2 = I_p*v2*I_p;
clear  u2 v2;

W = [ newu2'*conj(newv2), -newu2'*u1;...
    u1.'*conj(newu2), v1.'*u1 ];
W = (W-W.')/2; % -- here I am again manually anti-symmetrizate W for the
% same reason as before

Minv = inv( [ conj(newv2)/conj(newu2), -eye(N_row*N_col);...
    eye(N_row*N_col), -v1/u1 ] );

x_1 = min([x_i,x_f]);
y_1 = min([y_i,y_f]);
x_2 = max([x_i,x_f]);
y_2 = min([y_i,y_f]);
x_3 = max([x_i,x_f]);
y_3 = max([y_i,y_f]);

zeta = Minv(findex(x_1,y_1,N_col), findex(x_2,y_2,N_col)) ...
    + Minv(N_row*N_col + findex(x_1,y_1,N_col), findex(x_2,y_2,N_col)) ...
    - Minv(findex(x_1,y_1,N_col), N_row*N_col + findex(x_2,y_2,N_col)) ...
    - Minv(N_row*N_col + findex(x_1,y_1,N_col), N_row*N_col + findex(x_2,y_2,N_col))...
    + Minv(findex(x_2,y_2,N_col), findex(x_3,y_3,N_col)) ...
    + Minv(N_row*N_col + findex(x_2,y_2,N_col), N_row*N_col + findex(x_3,y_3,N_col)) ...
    - Minv(N_row*N_col + findex(x_2,y_2,N_col), findex(x_3,y_3,N_col)) ...
    - Minv(findex(x_2,y_2,N_col), N_row*N_col + findex(x_3,y_3,N_col));

op = (-1)^(N_row*N_col*(N_row*N_col+1)/2)*sqrt(abs(fdet(u1*lambda)))/fdet(u1*lambda)...
    *sqrt(abs(fdet(newu2*lambda)))/conj(fdet(newu2*lambda))*pfaffian_householder( W*lambda )...
    *zeta*(-1i);

save( strcat( "overlapY-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string( sp ),...
    ".mat" ),"op" );

% quit