function pnm = legendre(rlat, nmax)
%rlatlat:geocentric co-latitude
%Nmax: maximum expansion order
% pnm(1:nmax+1, 1:nmax+1):  associated legendre function
    pnm = zeros(nmax+1, nmax+1);
    z = cosd(rlat);

    p = plmbar( nmax, z);

    for n = 0:nmax
        for m = 0:n
            k = n*(n+1)/2 + m + 1;
            pnm(n+1, m+1) = p(k);

        end
    end
end
