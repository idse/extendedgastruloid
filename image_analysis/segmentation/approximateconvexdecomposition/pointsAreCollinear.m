function value = pointsAreCollinear(xy)
%check if a list of points stored as xy coordinates in an n x 2 array are
%collinear
%function is taken directly from John D'Errico's MATLAB answers response on
%6 January 2019

value = rank(xy(2:end,:) - xy(1,:)) == 1;
