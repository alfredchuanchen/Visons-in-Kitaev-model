load('parameters.mat');

% direc = "2e-braid/gauge_inv";
% direc = "2e-braid";

% direc = "4-site";

% direc = "3-site/upper";
direc = "3-site/lower";

% phase_list = zeros(L,1);

% --- The Berry phase from applying spin Z operators --- %

% for s = 1 : L
%     load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
%         direc,...
%         "/overlap-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
%         "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%         "-t-",string(t),...
%         "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%         "-z2-",string(z2),"-step-",string( s ),...
%         ".mat" ) );
%     phase_list(s) = op;
%     clear op;
% end
% 
% z = prod(phase_list);
% bp = imag( log( z/abs(z) ) );
% 
% % bp = 0;
% % 
% % for j = 1:L
% %     bp = bp + imag( log( phase_list(j)/abs(phase_list(j)) ) );
% % end
% 
% 
% save( strcat( "data/",direc,...
%     "/phase-2e",...
%     "-N1-",string(N_row),"-N2-",string(N_col),...
%     "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%     "-t-",string(t),...
%     "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%     "-z2-",string(z2),".dat" ), 'bp', '-ascii', '-double' );

% --- End --- %


% --- Berry phase from applying spin X operators --- %

% for s = 1 : L
%     load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
%         direc,...
%         "/overlapX-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
%         "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%         "-t-",string(t),...
%         "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%         "-z2-",string(z2),"-step-",string( s ),...
%         ".mat" ) );
% 
%     phase_list(s) = op;
%     clear op;
% end
% 
% z = prod(phase_list);
% bp = imag( log( z/abs(z) ) );
% 
% save( strcat( "data/",direc,...
%     "/phaseX-2e",...
%     "-N1-",string(N_row),"-N2-",string(N_col),...
%     "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%     "-t-",string(t),...
%     "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%     "-z2-",string(z2),".dat" ), 'bp', '-ascii', '-double' );

% --- End --- %

% --- Berry phase from applying X + Z --- %

% for s = 1 : L
%     load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
%         direc,...
%         "/overlap-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
%         "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%         "-t-",string(t),...
%         "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%         "-z2-",string(z2),"-step-",string( s ),...
%         ".mat" ) );
% 
%     phase_list(s) = op;
%     clear op;
%     
%     load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
%         direc,...
%         "/overlapX-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
%         "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%         "-t-",string(t),...
%         "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%         "-z2-",string(z2),"-step-",string( s ),...
%         ".mat" ) );
% 
%     phase_list(s) = phase_list(s) + op;
%     clear op;
% end
% 
% z = prod(phase_list);
% bp = imag( log( z/abs(z) ) );
% 
% save( strcat( "data/",direc,...
%     "/phaseZnX-2e",...
%     "-N1-",string(N_row),"-N2-",string(N_col),...
%     "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
%     "-t-",string(t),...
%     "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
%     "-z2-",string(z2),".dat" ), 'bp', '-ascii', '-double' );


% --- Berry phase for the upper / lower triangle --- %

L = 3;
phase_list = zeros(L,1);

load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
direc,"/overlapY-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
"-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
"-t-",string(t),"-Del_ep-",string(Del_ep),"-z1-",string(z1),...
"-z2-",string(z2),"-step-",string( 1 ),".mat" ) );

phase_list(1) = op;
clear op;

for s = 2 : L
    load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
        direc,...
        "/overlap-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
        "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
        "-t-",string(t),...
        "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
        "-z2-",string(z2),"-step-",string( s ),...
        ".mat" ) );

    phase_list(s) = op;
    clear op;
    
    load( strcat( "/Users/chuanchen/Dropbox/Chuan-Works/Toric code/fermion-flux/kitaev-honeycomb/data/",...
        direc,...
        "/overlapX-L-",string(L),"-N1-",string(N_row),"-N2-",string(N_col),...
        "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
        "-t-",string(t),...
        "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
        "-z2-",string(z2),"-step-",string( s ),...
        ".mat" ) );

    phase_list(s) = phase_list(s) + op;
    clear op;
end

if direc == "3-site/upper"

    z = prod(phase_list);
    bp = imag( log( z/abs(z) ) );
    name = "UT";

elseif direc == "3-site/lower"

    z = conj( prod(phase_list) );
    bp = imag( log( z/abs(z) ) );
    name = "LT";

end

save( strcat( "data/",direc,...
    "/phaseXYZ-",name,"-2e",...
    "-N1-",string(N_row),"-N2-",string(N_col),...
    "-J-",string(Jx),"-",string(Jy),"-",string(Jz),"-K-",string(K),...
    "-t-",string(t),...
    "-Del_ep-",string(Del_ep),"-z1-",string(z1),...
    "-z2-",string(z2),".dat" ), 'bp', '-ascii', '-double' );


% --- End --- %