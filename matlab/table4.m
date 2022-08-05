
nel = 128;
for p=8:16
    for k=0:p-1
        knot = [];
        for m=1:p-k % create multiple knots for reduced continuity
            knot = [knot, 0:nel];
        end
        for m=p-k:p % make sure start and end knot are repeated p+1 times
            knot = [knot, 0, nel];
        end
        knot = sort(knot);
        [w, x, rec, it] = getOptimalQuadPoints(knot, p);
        fprintf('%3d  ', rec);
    end
    fprintf('\n');
end
