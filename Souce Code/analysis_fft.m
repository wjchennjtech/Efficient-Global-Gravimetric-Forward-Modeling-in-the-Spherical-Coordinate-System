function [Cnm,Snm]=analysis_fft(data,nmax, dlat, dlon, nlat, nlon)


%Spherical Harmonic Analysis Using Fast Fourier Transform


Cnm=zeros(nmax+1,nmax+1);
Snm=zeros(nmax+1,nmax+1);

latmax=max(data(:,2));
latmin=min(data(:,2));
lonmin=min(data(:,1));
lonmax=max(data(:,1));
Lat=linspace(latmax,latmin,nlat);
Lon=linspace(lonmin,lonmax,nlon);

%Convert data to grid format
data_t=reshape(data(:,3),nlon,nlat);
data_t=data_t';

%Unified input of data based on the starting range of longitude
nlon2=nlon/2;
if lonmin < 0
    datat=zeros(size(data_t));
    datat(:,1:nlon2)=data_t(:,nlon2+1:nlon);
    datat(:,nlon2+1:nlon)=data_t(:,1:nlon2);
    data_t=datat;
end
data_t=fliplr(data_t);

%Angle to radians
dlat=dlat*pi/180;
dlon=dlon*pi/180;

for i=1:nlat

   %Obtain the Legendre function value corresponding to the latitude
    rlat=90-Lat(i);
    pnm = legendre(rlat, nmax);
    C=0;

    for j=1:nlon
        C(j)=complex(data_t(i,j)*dlon,0);
    end
    C2=fft(C, nlon);
    for n=0:nmax

        for m=0:n

            Cnm(n+1,m+1)=Cnm(n+1,m+1)+pnm(n+1,m+1)*real(C2(m+1))*sind(rlat)*dlat/4/pi;
            Snm(n+1,m+1)=Snm(n+1,m+1)+pnm(n+1,m+1)*imag(C2(m+1))*sind(rlat)*dlat/4/pi;

        end
    end
end


end

