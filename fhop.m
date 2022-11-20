function op = fhop(x,L)

if 0 < x && x <= L
    op = x;
elseif x > L
    op = mod(x,L);
else
    op = L+x;
end

end