function op = A_BdG(Nx,Ny,Jx,Jy,Jz,Del_ep,mate,z1,z2)

% we consider the type of bonds one by one, and we first construct the
% A matrix in related to the Majorana operators: H = I*A \gammga \gamma'

% --- x-bonds --- %

Ax = zeros(Nx*Ny);

for x = 1:Nx
    for y = 1:Ny-1
        Ax( findex(x,y,Ny), findex(fhop(x+1,Nx),y,Ny) ) = -Jx*(-1)^sum( mate(x,y+1:Ny) )...
            *z2^(x == Nx);
    end
%  right boundary no e parity term
    Ax(findex(x,Ny,Ny), findex(fhop(x+1,Nx),Ny,Ny)) = -Jx*z2^(x == Nx);
end

% --- y-bonds --- %

Ay = zeros(Nx*Ny);

for x = 1:Nx-1
    for y = 1:Ny-1
        Ay(findex(x,y,Ny), findex(x+1,y+1,Ny)) = -Jy*(-1)^(sum(mate(x,y+1:Ny)));
    end
    Ay( findex(x,Ny,Ny), findex(x+1,1,Ny)) = -Jy*(-1)^( sum(mate(x+1:Nx,:),'all') )...
        *z1;
end

for y = 1:Ny-1
    Ay(findex(Nx,y,Ny), findex(1,y+1,Ny)) = -Jy*(-1)^(sum(mate(Nx,y+1:Ny)))*z2;
end

Ay(findex(Nx,Ny,Ny), findex(1,1,Ny)) = -Jy*z1*z2;

% --- z-bonds --- %

Az = zeros(Nx*Ny);

for x = 1:Nx
    for y = 1:Ny-1
        Az( findex(x,y,Ny), findex(x,y+1,Ny) ) = -Jz;
    end
    Az( findex(x,Ny,Ny), findex(x,1,Ny) ) = -Jz*(-1)^( sum(mate(x:Nx,:),'all') )*z1;
end

A = Ax + Ay + Az + Del_ep/2*eye(Nx*Ny);

op = 2*[ 0, A;...
    -A.', 0 ];


end