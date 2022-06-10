
nel = 128;
for p=1:15
    for k=0:p-1
        knot = [];
        for m=1:p-k % create multiple knots for reduced continuity
            knot = [knot, 0:nel];
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
