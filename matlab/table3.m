
nel = 128;
alpha = 4/5;
for p=1:15
    for k=0:p-1
        knot = [];
        for m=1:p-k % create multiple knots for reduced continuity
            knot = [knot, alpha .^ [0:nel]];
        end
        for m=p-k:p % make sure start and end knot are repeated p+1 times
            knot = [knot, min(knot), max(knot)];
        end
        knot = sort(knot);
        [w, x, rec, it] = getOptimalQuadPoints(knot, p);
        if rec==1
            fprintf('%3d  ', it);
        else
            fprintf('  -  ');
        end
    end
    fprintf('\n');
end
