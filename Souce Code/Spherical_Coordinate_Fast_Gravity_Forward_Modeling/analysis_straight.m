function[result_cnm_snm]=analysis_straight(nmax, dlat, dlon, nlat, nlon,data)

%Directly perform spherical harmonic analysis

R=6371000;
rhoave=5500;

Cnm=zeros(nmax+1,nmax+1);
Snm=zeros(nmax+1,nmax+1);

CNM=zeros(nmax+1,nmax+1);
SNM=zeros(nmax+1,nmax+1);

for k=0:2  %Expand the quadratic term up to the third order
datah(:,1:2)=data(:,1:2);   
datah(:,3)=((data(:,4)./R).^(k+1))-((data(:,3)./R).^(k+1));
datah(:,3)=datah(:,3).*data(:,5);
[hcnm,hsnm]=analysis_fft(datah,nmax,dlat,dlon,nlat,nlon);

for l=0:nmax
CNM(l+1,:)=(3/(2*l+1))*(1/rhoave)*nchoosek(l+2,k)/(k+1).*(hcnm(l+1,:));
SNM(l+1,:)=(3/(2*l+1))*(1/rhoave)*nchoosek(l+2,k)/(k+1).*(hsnm(l+1,:));

end
Cnm=CNM+Cnm;
Snm=SNM+Snm;
end

nmnumber=0;
for i=1:nmax+1
    nmnumber=i+nmnumber;
end
result_cnm_snm=zeros(nmnumber,4);
row=1;

for n=0:nmax
    for m=0:n
        result_cnm_snm(row,1)=n;
        result_cnm_snm(row,2)=m;
        result_cnm_snm(row,3)=Cnm(n+1,m+1);
        result_cnm_snm(row,4)=Snm(n+1,m+1);
        row=1+row;
    end

end
end
