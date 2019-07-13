function x = tranpose_reshape(y)
[nRows,nCols,nPages]=size(y);

for j = 1:nCols
    for i = 1:nRows
        for k = 1:nPages
             x((j-1)*nPages*nRows+(i-1)*nPages+k)=y(i,j,k);
        end 
    end
end

% x = reshape(x_new,1,[]);