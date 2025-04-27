function p = plmbar(lmax, z)
%Calculate the value of Legendre function
    p = zeros((lmax+1)*(lmax+2)/2, 1);
    sqr = sqrt(1:2*lmax+1);

    f1 = zeros((lmax+1)*(lmax+2)/2, 1);
    f2 = zeros((lmax+1)*(lmax+2)/2, 1);

    phase = 1;
    scalef = 1.0e-280;

    k = 3;

    for l = 2:lmax
        k = k + 1;
        f1(k) = sqr(2*l-1) * sqr(2*l+1) / double(l);
        f2(k) = double(l-1) * sqr(2*l+1) / sqr(2*l-3) / double(l);
        for m = 1:l-2
            k = k+1;
            f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m);
            f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) / sqr(2*l-3) / sqr(l+m) / sqr(l-m);
        end
        k = k + 2;
    end

    u = sqrt((1.0-z)*(1.0+z));
    p(1) = 1.0;

    if lmax == 0
        return;
    end

    p(2) = sqr(3)*z;
    k = 2;

    for l = 2:lmax
        k = k + l;
        p(k) = f1(k)*z*p(k-l)-f2(k)*p(k-2*l+1);
    end

    pmm = sqr(2)*scalef;
    rescalem = 1.0/scalef;
    kstart = 1;

    for m = 1:lmax-1
        rescalem = rescalem * u;

        kstart = kstart+m+1;

        pmm = phase * pmm * sqr(2*m+1) / sqr(2*m);
        p(kstart) = pmm;

        k = kstart+m+1;
        p(k) = z * sqr(2*m+3) * pmm;

        for l = m+2:lmax
            k = k+l;
            p(k) = z*f1(k)*p(k-l)-f2(k)*p(k-2*l+1);
            p(k-2*l+1) = p(k-2*l+1) * rescalem;
        end

        p(k) = p(k) * rescalem;
        p(k-lmax) = p(k-lmax) * rescalem;
    end

    rescalem = rescalem * u;
    kstart = kstart+m+1;
    p(kstart) = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax) * rescalem;

end
