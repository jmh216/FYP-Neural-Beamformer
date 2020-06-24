function[axh] = fcn_20171103_01_subplots(nRow,nCol)
% quick function to make a subplot grid where handles has same dimensions
% as grid

id = 0;
for icol = 1:nCol
    for irow=1:nRow
        id = id+1;
        axh(id) = subplot(nRow,nCol,(irow-1)*nCol + icol);
    end
end
axh = reshape(axh,nRow,nCol);