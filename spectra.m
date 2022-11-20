% ---import the parameterization of the model--- %
load( 'parameters.mat' );

% tic

sp = 0;

% Tune the e particle parity:

% ---Move one e particle away horizontally--- %

% if sp == 0
% else
%     mate( fix(N_row/2)+1, 1 ) = 1;
% % ---C = 1--- %
%     mate( path(sp, 1), path(sp, 2) ) = 1;
% %  ---C = 2--- %
% %     mate( path(sp, 1), path(sp, 2) ) = -1;
% % ---end--- %
% end


% ---the 4-site / 3-site model--- %
% we fix one e-particle at the N_row/2 row and first column

if sp == 0
else
    mate( fix(N_row/2)+1, 1 ) = 1;
% ---C = 1--- %
    mate( path(sp, 1), path(sp, 2) ) = 1;
% ---C = 2--- %
%     mate( path(sp, 1), path(sp, 2) ) = -1;
end


% ---braiding of two e particles (closed loop)--- %
% we fix one of the e particles at [ fix(N_row/2), fix(N_col/2) ]
% 
% if sp == 0
% else
%     mate( fix(N_row/2), fix(N_col/2) ) = 1;
%     mate( path(sp,1), path(sp,2) ) = 1;
% end

% ---braiding in a gauge invariant way--- %

% mate( path(sp,1), path(sp,2) ) = 1;
% % 
% % mate( path(sp,3), path(sp,4) ) = 1;
% % ---For C = 2, the local potential are different for the 2 e particles--- %
% mate( path(sp,3), path(sp,4) ) = -1;
% % ---end--- %


% ---start the diagonalization--- %

h_BdG = H_BdG(N_row,N_col,Jx,Jy,Jz,K,t,Del_ep,mate,z1,z2);

[ eigenvec, eigenval ] = eig( h_BdG );

[ eigenval, ind ] = sort( diag(eigenval), 'descend' );
eigenvec = eigenvec( :, ind );

% We ultimately want the energies sligned in the form { E_1, E_2, ... ,
% -E_1, -E_2, ... , -E_N }, with E_1 < E_2 < ... < E_N. We need to sort
% the positive energies in an ascending order:

[ poseigenval, ind1 ] = sort( eigenval( 1 : N_row*N_col ), 'ascend' );
eigenvec( :, 1 : N_row*N_col) = eigenvec( :, ind1 );

eigenval( 1 : N_row*N_col ) = poseigenval;

% --done-- %


% ---define the new eigenvector matrix in the form of [ u, v^*; v, u^* ]:
% this is because the eigenvectors calculated numerically may differ by
% some phase factor like "-1"

neweigenvec = zeros( 2*N_row*N_col );
neweigenvec( 1:2*N_row*N_col, 1:N_row*N_col )...
    = eigenvec( 1:2*N_row*N_col, 1:N_row*N_col );
neweigenvec( 1:N_row*N_col, N_row*N_col+1 : 2*N_row*N_col )...
    = conj( eigenvec( N_row*N_col+1 : 2*N_row*N_col, 1:N_row*N_col ) );
neweigenvec( N_row*N_col+1 : 2*N_row*N_col, N_row*N_col+1 : 2*N_row*N_col )...
    = conj( eigenvec( 1 : N_row*N_col, 1:N_row*N_col ) );

% du = det( eigenvec( 1:N_row*N_col, 1:N_row*N_col ) );
% dv = det( eigenvec( N_row*N_col+1:2*N_row*N_col, 1:N_row*N_col ) );

% ---export the data--- %

% save( strcat( "spectra-N1-",string(N_row),"-N2-",string(N_col),...
%     "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%     "-t-",string(t),...
%     "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%     "-z2-",string(z2),"-step-",string(step),...
%     ".mat" ),"neweigenvec", "mate" );

save( strcat( "spectra-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string(sp),...
    ".dat" ), "eigenval", '-ascii', '-double' );

% lowvec = neweigenvec(:,N_col*N_row); % the lowest energy eigenvector of H_BdG
% 
% save( strcat( "spectra-N1-",string(N_row),"-N2-",string(N_col),...
%     "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%     "-t-",string(t),...
%     "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%     "-z2-",string(z2),"-step-",string(step),...
%     ".mat" ), "lowvec" );


% ---calculate the gound state parity--- %

P_mat = zeros(2*N_row*N_col);
for i = 1:N_row*N_col
    P_mat(i, 2*i-1) = 1;
    P_mat(N_row*N_col+i, 2*i) = 1;
end
% permutate the (u_1, ..., u_N, v_1, ..., v_N) into 
% (u_1, v_1, ... u_N, v_N)
% note here we have E_1 < E_2 < E_3 < ... < E_N

OD = 1/sqrt(2)*[ eye(N_row*N_col), eye(N_row*N_col);...
    -1i*eye(N_row*N_col), 1i*eye(N_row*N_col) ]*neweigenvec*P_mat*...
    kron( eye(N_row*N_col), [ 0, 1; 1, 0 ] );
B = 1i*OD*kron( eye(N_row*N_col), [1,0;0,-1] )*OD';


% use different packages to calculate the Pfaffian:
% parity = pfaffian_LTL( B )*sqrt(2)*cos(N_row*N_col*pi/2-pi/4);
parity = pfaffian_householder( B )*sqrt(2)*cos(N_row*N_col*pi/2-pi/4);
parity = round(parity);
% the final cosine factor is due to the difference in defining the fermions.

E_GS = 1/2*( sum(eigenval(N_row*N_col+1:2*N_row*N_col) ) );

% the average fermion number at plaquette (N_row/2+1, N_col/2)
n_ave = sum( neweigenvec( N_row*N_col+findex( N_row/2+1, N_col/2, N_col ), 1:N_row*N_col ).*...
    conj(neweigenvec( N_row*N_col+findex( N_row/2+1, N_col/2, N_col ), 1:N_row*N_col )) );

save( strcat( "parity-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string(sp),...
    ".dat" ), "parity", "E_GS", "n_ave", '-ascii', '-double' );

% Now, depends on the ground state parity, if the GS parity is odd, the state
% that we are interested will be the first excited state of that Hamiltonian, \Psi,
% and we need to define a new Hamiltonian such that \Psi is the ground
% state of this newly defined Hamiltonian.
% And in this new Hamiltonian, the newly defined [ u, v^*; v, u^* ] are
% related to those in the original Hamiltonian.

u1 = neweigenvec( 1 : N_row*N_col, 1 );
v1 = neweigenvec( 1+N_row*N_col : 2*N_row*N_col, 1 );

u2 = neweigenvec( 1 : N_row*N_col, 2 );
v2 = neweigenvec( 1+N_row*N_col : 2*N_row*N_col, 2 );

poseigenvec = neweigenvec(:, 1:N_row*N_col);

if parity == -1
    poseigenvec( 1:N_row*N_col, 1 ) = conj(v1);
    poseigenvec( 1+N_row*N_col:2*N_row*N_col, 1 ) = conj(u1);
else
end

% ---export all the single particle modes' eigenvectors--- %

save( strcat( "spectra-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),"-step-",string(sp),...
    ".mat" ),"poseigenvec", "mate");


% ---export the two lowest energy single particle modes--- %

% save( strcat( "LE-mode-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
%     "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%     "-t-",string(t),...
%     "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%     "-z2-",string(z2),"-step-",string(sp),...
%     ".mat" ), "u1", "v1", "u2", "v2", "mate");

% toc

% quit