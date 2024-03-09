function r = drawmnrnd(n,p)
%n is an integer specifying the number of draws we want to perform from a
%multinomial distribution
%p is a vector of probabilities
%r is an nx1 vector with the resulting random numbers


p = p/sum(p); %normalize p (e.g. in case of rounding errors in specifying probabilities)
nv = length(p);
cp = cumsum(p);

r = NaN(n,1);
q = rand(n,1);
for ii = 1:n
    for jj = 1:nv
        if q(ii) <= cp(jj)
            r(ii) = jj;
            break
        end
    end
end


end



