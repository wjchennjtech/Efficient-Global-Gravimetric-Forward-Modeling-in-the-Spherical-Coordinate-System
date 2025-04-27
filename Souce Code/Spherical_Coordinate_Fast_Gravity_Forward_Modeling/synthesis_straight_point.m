function[gravity_potential,gravity_anomal,gravity_gradient]=synthesis_straight_point(cnm_snm,nmax,Lat,Lon,height)

%Point data spherical harmonic synthesis
R=6371000;
GM=3986004.415*1e+8;
row=size(Lat,1);

gravity_potential=zeros(row,1);
gravity_anomal=zeros(row,1);
gravity_gradient=zeros(row,1);
rp=height+R;

Cnm=zeros(nmax+1,nmax+1);
Snm=zeros(nmax+1,nmax+1);

nmnumber=0;
for i=1:nmax+1
    nmnumber=i+nmnumber;
end

for i=1: nmnumber

    Cnm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,3);
    Snm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,4);

end

for i=1:row
        rlat=90-Lat(i);
        pnm = legendre(rlat, nmax);    
for n=0:nmax

    for m=0:n
        coef=GM/R;

        potential=(R/rp).^(n+1);
        anomal=(n+1).*potential./rp;
        gradient=(n+1).*(n+2).*potential./(rp*rp);

        gravity_potential(i)=coef.*potential.*pnm(n+1,m+1)*Cnm(n+1,m+1).*cosd(m*Lon(i))...
            +coef.*potential.*pnm(n+1,m+1)*Snm(n+1,m+1).*sind(m*Lon(i))...
            +gravity_potential(i);
        gravity_anomal(i)=coef.*anomal.*pnm(n+1,m+1)*Cnm(n+1,m+1).*cosd(m*Lon(i))...
            +coef.*anomal.*pnm(n+1,m+1)*Snm(n+1,m+1).*sind(m*Lon(i))...
            +gravity_anomal(i);
        gravity_gradient(i)=coef.*gradient.*pnm(n+1,m+1)*Cnm(n+1,m+1).*cosd(m*Lon(i))...
            +coef.*gradient.*pnm(n+1,m+1)*Snm(n+1,m+1).*sind(m*Lon(i))...
            +gravity_gradient(i);

    end

end
end
