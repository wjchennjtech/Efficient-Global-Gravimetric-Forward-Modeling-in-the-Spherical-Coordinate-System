function [gravity_potential,gravity_anomal,gravity_gradient]=synthesis_straight(cnm_snm,nmax,latmax,latmin,lonmax,lonmin,dlat,dlon,height)

R=6371000;
GM=3986004.415*1e+8;

%The angle of each latitude and longitude band
nlat=180/dlat;
nlon=360/dlon;
Lat=linspace(latmax,latmin,nlat);
Lon=linspace(lonmin,lonmax,nlon);
rp=height+R;%Computation height

gravity_potential=zeros(nlat,nlon);
gravity_anomal=zeros(nlat,nlon);
gravity_gradient=zeros(nlat,nlon);

Cnm=zeros(nmax+1,nmax+1);
Snm=zeros(nmax+1,nmax+1);

%Extraction coefficient
nmnumber=0;
for i=1:nmax+1
    nmnumber=i+nmnumber;
end

for i=1: nmnumber

    Cnm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,3);
    Snm(cnm_snm(i,1)+1,cnm_snm(i,2)+1)=cnm_snm(i,4);

end

for i=1:nlat

  %Obtain the Legendre function value corresponding to the latitude
    rlat = 90 - Lat(i);
    pnm = legendre(rlat, nmax);

    CS_p=zeros(nmax+1);
    CS_a=zeros(nmax+1);
    CS_g=zeros(nmax+1);
    for m=0:nmax
        alpham_p=0;beltam_p=0;
        alpham_a=0;beltam_a=0;
        alpham_g=0;beltam_g=0;
        for n=m:nmax
            potential=(R/rp).^(n+1);
            anomal=(n+1).*potential./rp;
            gradient=(n+1).*(n+2).*potential./(rp*rp);

            alpham_p=alpham_p+Cnm(n+1,m+1)*pnm(n+1,m+1)*potential;
            beltam_p=beltam_p+Snm(n+1,m+1)*pnm(n+1,m+1)*potential;
            alpham_a=alpham_a+Cnm(n+1,m+1)*pnm(n+1,m+1)*anomal;
            beltam_a=beltam_a+Snm(n+1,m+1)*pnm(n+1,m+1)*anomal;
            alpham_g=alpham_g+Cnm(n+1,m+1)*pnm(n+1,m+1)*gradient;
            beltam_g=beltam_g+Snm(n+1,m+1)*pnm(n+1,m+1)*gradient;

        end
        cm_p=alpham_p/2;sm_p=beltam_p/2;
        cm_a=alpham_a/2;sm_a=beltam_a/2;
        cm_g=alpham_g/2;sm_g=beltam_g/2;
        CS_p(m+1)=complex(cm_p,-sm_p);
        CS_a(m+1)=complex(cm_a,-sm_a);
        CS_g(m+1)=complex(cm_g,-sm_g);
    end
    CS_p2=fft(CS_p, nlon);
    CS_a2=fft(CS_a, nlon);
    CS_g2=fft(CS_g, nlon);
    gravity_potential(i,1:nlon)=2*real(CS_p2(1:nlon));
    gravity_anomal(i,1:nlon)=2*real(CS_a2(1:nlon));
    gravity_gradient(i,1:nlon)=2*real(CS_g2(1:nlon));
end

%Adjust output data based on the starting range of longitude
nlon2=nlon/2;
if lonmin < 0
    fp=zeros(size(gravity_potential));
    fp(:,1:nlon2)=gravity_potential(:,nlon2+1:nlon);
    fp(:,nlon2+1:nlon)=gravity_potential(:,1:nlon2);
    gravity_potential=fp;
    fa=zeros(size(gravity_anomal));
    fa(:,1:nlon2)=gravity_anomal(:,nlon2+1:nlon);
    fa(:,nlon2+1:nlon)=gravity_anomal(:,1:nlon2);
    gravity_anomal=fa;
    fg=zeros(size(gravity_gradient));
    fg(:,1:nlon2)=gravity_gradient(:,nlon2+1:nlon);
    fg(:,nlon2+1:nlon)=gravity_gradient(:,1:nlon2);
    gravity_gradient=fg;
end

coef=GM/R;
gravity_potential=gravity_potential.*coef;
gravity_anomal=gravity_anomal.*coef;
gravity_gradient=gravity_gradient.*coef;
gravity_potential=fliplr(gravity_potential);
gravity_anomal=fliplr(gravity_anomal);
gravity_gradient=fliplr(gravity_gradient);

end