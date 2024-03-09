function val = doIntersect(p1, q1, p2, q2)
%check if the line segments p1,q1 and p2,q2 intersect
%return true if so, false if not

o1 = orientation(p1, q1, p2);
o2 = orientation(p1, q1, q2);
o3 = orientation(p2, q2, p1);
o4 = orientation(p2, q2, q1);

if (o1 ~= o2) && (o3 ~= o4)
    val = true;
elseif (o1 == 0) && onSegment(p1, p2, q1)
    val = true;
elseif (o2 == 0) && onSegment(p1, q2, q1)
   val = true;
elseif (o3 == 0) && onSegment(p2, p1, q2)
    val = true;
elseif (o4 == 0) && onSegment(p2, q1, q2)
    val = true;
else
    val = false;    
end

end