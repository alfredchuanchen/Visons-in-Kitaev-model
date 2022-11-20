% In this script, we are calculating the overlap between two states with a
% superposition of different amount of hole pairs on top of the fully
% occupied state
% This is mainly used for certain regime of the parameters, e.g., Kitaev
% model with Jx = 1, Jy = -1, Jz = 1

load('parameters.mat');

sp = 2;

% lambda = 2.85 - 1.85*tanh(abs(Del_ep)/1.1); % for FM case with small K
% lambda = 2.21 - 1.21*tanh(abs(Del_ep)/1.3); % for FM case

% lambda = 1.61 - 0.61*tanh(abs(Del_ep)/2); % for AFM case
lambda = 1.61;
% this is a parameter used to increase the determinant and
% the pfaffian to to avoid the numerical 0.


% ---when moving only one particle, the row and column coordinate along the path---%

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


% this is the "bra" state (the state after acting the Z operaotr)

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

% this is the ket state (the state before acting the Z operator)

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

if x_f == x_i

    % in this case the hopping of "e" is along the horizontal direction, so
    % does not cross the braanch cut of fermion (as long as not at the
    % boundary)
    
    W = [ v2'*conj(u2), -v2'*v1;...
        v1.'*conj(v2), u1.'*v1 ];
    W = (W-W.')/2; % -- here I am manually anti-symmetrizate M to same
    % time for the numerical check of whether M is skew-symmetric or not in
    % the pfaffian_LTL function

    op = (-1)^(N_row*N_col*(N_row*N_col+1)/2)*sqrt( abs(fdet(v1*lambda)) )/fdet(v1*lambda)...
        *sqrt(abs(fdet(v2*lambda)))/conj(fdet(v2*lambda))*pfaffian_householder( W*lambda );
    % one should be careful about the power on top of '-1' in the formula
    % we used, should be the linear size of the matrix u, which is N_row*N_col,
    % i.e., the system size!
else
    I_p = eye(N_row*N_col);
    for y = 1:y_i-1
        I_p( findex( max([ x_i, x_f ]), y, N_col ), findex( max([ x_i, x_f ]), y, N_col ) ) = -1;
    end
    newu2 = I_p*u2*I_p;
    newv2 = I_p*v2*I_p;
    clear  u2 v2;

    W = [ newv2'*conj(newu2), -newv2'*v1;...
        v1.'*conj(newv2), u1.'*v1 ];
    W = (W-W.')/2; % -- here I am again manually anti-symmetrizate M for the
    % same reason as before
    op = (-1)^(N_row*N_col*(N_row*N_col+1)/2)*sqrt(abs(fdet(v1*lambda)))/fdet(v1*lambda)...
        *sqrt(abs(fdet(newv2*lambda)))/conj(fdet(newv2*lambda))*pfaffian_householder( W*lambda )...
        *(-1)^(y_i-1);
    
end

save( strcat( "overlap-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string( sp ),...
    ".mat" ),"op" );

% quit