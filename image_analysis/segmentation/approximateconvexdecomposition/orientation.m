function val = orientation(p, q, r)
%return the orientation (clockwise, counterclockwise, or colinear) of an
%ordered triplet of points p, q, and r. Each input point should be a
%two-element vector with x and y coordinates of that point in the cartesian
%plane

test = (q(2) - p(2))*(r(1) - q(1)) - (q(1) - p(1))*(r(2) - q(2));

if test > 0
    %clockwise
    val = 1;
elseif test < 0
    %counterclockwise
    val = 2;
else
    %colinear
    val = 0;
end

end