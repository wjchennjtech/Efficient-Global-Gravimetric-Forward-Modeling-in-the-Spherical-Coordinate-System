function[Cnm,Snm]=analysis_slice(nmax, dlat, dlon, nlat, nlon, xyz, nz)

%Layered spherical harmonic analysis

R=6371000;
rhoave=5500;

Cnm=zeros(nmax+1,nmax+1,nz);
Snm=zeros(nmax+1,nmax+1,nz);

[CSnm,z] = divide_layer(xyz,nz,dlat,dlon,nlat,nlon);

%Perform spherical harmonic analysis on each layer
for i=1:nz

    data=CSnm(:,:,i);

    [cnm,snm]=analysis_fft(data,nmax,dlat,dlon,nlat,nlon);
    CNM=zeros(nmax+1,nmax+1);
    SNM=zeros(nmax+1,nmax+1);
    
    for n=0:nmax
        z1=(z(i+1)/R).^(n+3);
        z2=(z(i)/R).^(n+3);
        zdiff=z1-z2;
        coef=3/rhoave/(2*n+1)/(n+3);
        CNM(n+1,:)=cnm(n+1,:)*zdiff*coef;
        SNM(n+1,:)=snm(n+1,:)*zdiff*coef;
    end

    Cnm(:,:,i)=CNM(:,:);
    Snm(:,:,i)=SNM(:,:);

end

end
