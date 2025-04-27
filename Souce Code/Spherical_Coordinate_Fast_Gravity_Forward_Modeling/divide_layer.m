function [CSNM,z] = divide_layer(data,nz,dlat,dlon,nlat,nlon)
% divide_layer

R=6371*1e+3;

lnum=length(data);

zmax=max(data(:,4));
zmin=min(data(:,3));
latmin=min(data(:,2));latmax=max(data(:,2));
lonmin=min(data(:,1));lonmax=max(data(:,2));

CSnm=zeros(nlat,nlon,nz);

z=linspace(zmin,zmax,nz+1);
z=z+R;
dz=abs(z(1)-z(2));

%Store in layers according to needs for easier future access
for i=1:size(data, 1)
    l1=floor(ceil((data(i,3)-zmin)/(dz/2))/2); 
    l2=floor(ceil((data(i,4)-zmin)/(dz/2))/2);
    l=(l1:l2);
    l(find(l==0))=[];
    lon = ceil((data(i,1)-lonmin)/dlon)+1; 
    lat = ceil(latmax-(data(i,2))/dlat)+1;
    CSnm(lat,lon,l) = data(i,5);    
end

CSNM=zeros(nlat*nlon,3);
for n=1:nz
    for i=1:nlat
        for j=1:nlon
            CSNM((i-1)*nlon+j,1,n)=(j-1)*dlon+lonmin;
            CSNM((i-1)*nlon+j,2,n)=latmin+(i-1)*dlat;
            CSNM((i-1)*nlon+j,3,n)=CSnm(i,j,n);
        end
    end
end
