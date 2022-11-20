% ---parameter of the model--- %

N_row = 24; % number of rows, i.e., along the verticle direction
N_col = 24; % number of columns, i.e., along the horizontal direction

Jx = -1;
Jy = -1;
Jz = -1;

Del_ep = 0;
% epsilon onsite energy: \Delta_epsilon c^\dagger c

K = 0.01;
% Haldane mass coming from the external field, K = h^3/J^2

t = 0.;

z1 = -1; % boundary condition along the horizontal direction
z2 = -1; % BC along the verticle direction

mate = zeros(N_row,N_col);


% ---the 4-site model's path--- %

% L = 4;
% 
% path = [ fix(N_row/2)+1, fix(N_col/2);
%     fix(N_row/2)+1, fix(N_col/2)+1;
%     fix(N_row/2), fix(N_col/2)+1;
%     fix(N_row/2), fix(N_col/2) ];


% --- Upper/lower triangle path --- %

L = 3;

pname = "upper";
% pname = "lower";

if pname == "upper"

    path = [ fix(N_row/2)+1, fix(N_col/2);
    fix(N_row/2), fix(N_col/2) + 1;
    fix(N_row/2), fix(N_col/2) ];

elseif pname == "lower"

    path = [ fix(N_row/2)+1, fix(N_col/2);
    fix(N_row/2), fix(N_col/2)+1;
    fix(N_row/2)+1, fix(N_col/2)+1 ];

end


% ---Movement of one e particle away horizontally--- %

% in this case, we place 1 of the e particle at the first vertice of a row,
% and move the other one from the site next to it till the boundary

% L = N_col-1;
% path = [ ( fix(N_row/2)+1 )*ones(L,1), linspace(2,N_col,L).' ];



% ---braiding two e particles--- %
% the position of the fixed e particle would be [ fix(N_row/2), fix(N_col/2) ]

% ---the contour with even number of enclosed plaquette--- %
% 
% r = fix(N_row/4);
% 
% path1 = [ linspace(fix(N_row/2)-r, fix(N_row/2)+r-1, 2*r).',...
%     ( fix(N_col/2)-r )*ones(2*r,1) ];
% 
% path2 = [ (fix(N_row/2)+r )*ones(2*r,1),...
%     linspace(fix(N_col/2)-r, fix(N_col/2)+r-1, 2*r).' ];
% 
% path3 = [ linspace(fix(N_row/2)+r, fix(N_row/2)-r+1, 2*r).',...
%     ( fix(N_col/2)+r )*ones(2*r,1) ];
% 
% path4 = [ (fix(N_row/2)-r )*ones(2*r,1),...
%     linspace(fix(N_col/2)+r, fix(N_col/2)-r+1, 2*r).' ];
% 
% path = [ path1; path2; path3; path4 ];
% 
% L = size(path,1);


% ---contour with odd number of plaquette--- %

% r = fix(N_row/4);
% 
% path1 = [ linspace(fix(N_row/2)-r, fix(N_row/2)+r, 2*r+1).',...
%     ( fix(N_col/2)-r )*ones(2*r+1,1) ];
% 
% path2 = [ (fix(N_row/2)+r+1 )*ones(2*r+1,1),...
%     linspace(fix(N_col/2)-r, fix(N_col/2)+r, 2*r+1).' ];
% 
% path3 = [ linspace(fix(N_row/2)+r+1, fix(N_row/2)-r+1, 2*r+1).',...
%     ( fix(N_col/2)+r+1 )*ones(2*r+1,1) ];
% 
% path4 = [ (fix(N_row/2)-r )*ones(2*r+1,1),...
%     linspace(fix(N_col/2)+r+1, fix(N_col/2)-r+1, 2*r+1).' ];
% 
% path = [ path1; path2; path3; path4 ];
% 
% L = size(path,1);



% ---"gauge invariant" path for braiding--- %

% % the four end points of the path
% %
% x1 = fix(N_row/2);
% y1 = fix(N_col/4) - 0;
% % here we have changed the initial position of one particle in horizontal
% % direction by 1.
% 
% x2 = fix(N_row/2);
% y2 = fix(N_col/2);
% 
% x3 = fix(N_row/2);
% y3 = fix(N_col*3/4);
% 
% x4 = fix(N_row/4);
% y4 = fix(N_col/2);

% there will be 5 intermediate paths for doing the braiding

% nada = linspace(x2,x4+1,x2-x4).';
% path1 = [ x1*ones(size(nada)), y1*ones(size(nada)), nada,...
%     y2*ones(size(nada)) ];
% 
% nada = linspace(y1,y3-1,y3-y1).';
% path2 = [ x1*ones(size(nada)), nada, x4*ones(size(nada)), y4*ones(size(nada)) ];
% 
% nada = linspace(x4,x2-1,x2-x4).';
% path3 = [ x3*ones(size(nada)), y3*ones(size(nada)), nada, y2*ones(size(nada)) ];
% 
% nada = linspace(y2,y1+1,y2-y1).';
% path4 = [ x3*ones(size(nada)), y3*ones(size(nada)), x1*ones(size(nada)), nada ];
% 
% nada = linspace(y3,y2+1,y3-y2).';
% path5 = [ x2*ones(size(nada)), nada, x1*ones(size(nada)), y1*ones(size(nada)) ];
% 
% % path = [ path1; path2; path3; path4; path5 ];
% 
% %---for double exchange with the C = 2 case---% 
% 
% nada = linspace(x2,x4+1,x2-x4).';
% path6 = [ nada, y2*ones(size(nada)), x1*ones(size(nada)), y1*ones(size(nada)) ];
% 
% nada = linspace(y1,y3-1,y3-y1).';
% path7 = [ x4*ones(size(nada)), y4*ones(size(nada)), x1*ones(size(nada)), nada ];
% 
% nada = linspace(x4,x2-1,x2-x4).';
% path8 = [ nada, y2*ones(size(nada)), x3*ones(size(nada)), y3*ones(size(nada)) ];
% 
% nada = linspace(y2,y1+1,y2-y1).';
% path9 = [ x1*ones(size(nada)), nada, x3*ones(size(nada)), y3*ones(size(nada)) ];
% 
% nada = linspace(y3,y2+1,y3-y2).';
% path10 = [ x1*ones(size(nada)), y1*ones(size(nada)), x3*ones(size(nada)), nada ];
% 
% path = [ path1; path2; path3; path4; path5; path6; path7; path8; path9; path10 ];
% 
% L = size(path,1);


% ---export parameters--- %

disp(strcat("-> The number of total steps is L = ",string(L)));

save( 'parameters.mat', 'N_row', 'N_col', 'Jx', 'Jy', 'Jz','K',...
    't', 'Del_ep','z1', 'z2', 'path', 'mate', 'L' ); 


% --- without specifying the K --- %

% save( 'parameters.mat', 'N_row', 'N_col', 'Jx', 'Jy', 'Jz', ...
%     't', 'Del_ep','z1', 'z2', 'path', 'mate', 'L' );


% --- without specifying the Del_ep and K --- %

% save( 'parameters.mat', 'N_row', 'N_col', 'Jx', 'Jy', 'Jz',...
%     't', 'z1', 'z2', 'path', 'mate', 'L' ); 


