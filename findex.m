function nop = findex(x,y,Ny)
% here the x is the row index (from top to bottom) and y is the column
% index (from left to right), Ny is the number of columns

nop = y + (x-1)*Ny;

end