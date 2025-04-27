function[data,nlayers]=transform_data(table_data1,table_data2)

data=[];
nlayers=size(table_data1,2)-3;
for i=1:nlayers
    dataup=table_data1(:,[1,2,i+2]);
    datalow=table_data1(:,[1,2,i+3]);
    deltarho=table_data2(:,[1,2,i+2]);
    [data0]=blend_data(datalow,dataup,deltarho);
    data=[data;data0];
end
end


function[data_b]=blend_data(datalow,dataup,deltarho)
data_b(:,[1,2])=dataup(:,[1,2]);
data_b(:,3)=datalow(:,3)*1e+3;
data_b(:,4)=dataup(:,3)*1e+3;
data_b(:,5)=deltarho(:,3)*1e+3;
end