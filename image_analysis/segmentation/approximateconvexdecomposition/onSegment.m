function val = onSegment(p, q, r)
%for three colinear input points p, q, r, check if q lies between p and r

if (q(1) <= max(p(1), r(1))) && (q(1) >= min(p(1), r(1))) && (q(2) <= max(p(2), r(2))) && (q(2) >= min(p(2),r(2)))
    val = true;
else
    val = false;
end


end