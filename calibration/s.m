cl=0;
if(cl==1)
  clear all;
  close all;
end
part=2;
aks2_=load('ins-mems-aks2.log','-ascii');
aks1_=load('ins-mems-aks1.log','-ascii');
dus2_= load('ins-mems-dus2.log','-ascii');

a=aks1_(:,5:7);
b=dus2_(:,8:10);
c=aks2_(:,5:7);

if(part==1)
  figure;
  plot(a(:,1),'r');
  hold on;
  plot(a(:,2),'g');
  plot(a(:,3),'b');

  figure;
  plot(b(:,1),'r');
  hold on;
  plot(b(:,2),'g');
  plot(b(:,3),'b');


  figure;
  plot(c(:,1),'r');
  hold on;
  plot(c(:,2),'g');
  plot(c(:,3),'b');
end
if(part==2)
  n=6;
  t1=[1,700,1200,1700,2200,2700];
  t2=[500,1100,1500,2100,2500,3200];
  for i=1:n
    s1(i)=sum(a(t1(i):t2(i),1))/(t2(i)-t1(i)+1);
    s2(i)=sum(a(t1(i):t2(i),2))/(t2(i)-t1(i)+1);
    s3(i)=sum(a(t1(i):t2(i),3))/(t2(i)-t1(i)+1);
  end
  s1
  s2
  s3
  
  ba(1)=1/2 *(s1(1)+s1(2));  
  ba(2)=1/2 *(s2(1)+s2(2));  
  ba(3)=1/2 *(s3(1)+s3(2));
  
  Sa(1,:)=1/2 /9.8 * [(s1(6)-s1(5)); (s1(4)-s1(3)); (s1(1)-s1(2))];
  Sa(2,:)=1/2 /9.8 * [(s2(6)-s2(5)); (s2(4)-s2(3)); (s2(1)-s2(2))];
  Sa(3,:)=1/2 /9.8 * [(s3(6)-s3(5)); (s3(4)-s3(3)); (s3(1)-s3(2))];
  
  Sa

  Sf=inv(-Sa);
  bf=-Sf*ba';
  
  Sf
  bf
  
  fz=Sf*a'+bf;
  fz=fz';
  
  figure;
  plot(fz);
end%if
if(part==3)
  n=6;
  t1=[50,1000,2050,2950,4350,5100];
  t2=[600,1450,2600,3400,4750,5500];
  b_w=[sum(b(3500:4000,1));sum(b(3500:4000,2));sum(b(3500:4000,3))]/501;
  b(:,1)=b(:,1)-b_w(1);
  b(:,2)=b(:,2)-b_w(2);
  b(:,3)=b(:,3)-b_w(3);
  for i=1:n
    dt=dus2_(t1(i),12)*60 + dus2_(t1(i),13) + dus2_(t1(i),14)*0.001 - (dus2_(t2(i),12)*60 + dus2_(t2(i),13) + dus2_(t2(i),14)*0.001);
    s1(i)=dt*sum(b(t1(i):t2(i),1))/(t2(i)-t1(i)+1);
    s2(i)=dt*sum(b(t1(i):t2(i),2))/(t2(i)-t1(i)+1);
    s3(i)=dt*sum(b(t1(i):t2(i),3))/(t2(i)-t1(i)+1);
  end
  S_w(1,:)=1/2 /180 * [(s1(6)-s1(5)); (s1(3)-s1(4)); (s1(1)-s1(2))];
  S_w(2,:)=1/2 /180 * [(s2(6)-s2(5)); (s2(3)-s2(4)); (s2(1)-s2(2))];
  S_w(3,:)=1/2 /180 * [(s3(6)-s3(5)); (s3(3)-s3(4)); (s3(1)-s3(2))];
  
  S_om=-inv(S_w);

  S_om
  
  fi=S_om*b';
  
  figure;
  plot(fi(1,:),'r');
  hold on;
  plot(fi(2,:),'g');
  plot(fi(3,:),'b');
end
if (part==4) 
  n=20;
  t1=[50,400,800,1400,1900,2300,2800,3350,3750,4200,4550,5050,5500,5950,6300,6700,7100,7550,8000,8350];%íà÷àëà
  t2=[250,700,1250,1750,2200,2600,3000,3600,4050,4450,4900,5350,5800,6200,6550,6950,7400,7900,8200,8498];%êîíöû
  for i=1:n
    x(i)=sum(c(t1(i):t2(i),1))/(t2(i)-t1(i));
    y(i)=sum(c(t1(i):t2(i),2))/(t2(i)-t1(i));
    z(i)=sum(c(t1(i):t2(i),3))/(t2(i)-t1(i));
  end%for
  
  x2=x.^2;
  y2=y.^2;
  z2=z.^2;
  xy=x.*y.*2;
  zx=z.*x.*2;
  yz=y.*z.*2;
  H=[x2;y2;z2;xy;zx;yz;x;y;z];
  MM=-inv(H*(H'))*H*ones(n,1);
  
  mt=MM(7:9);
  Mt(1,1)=MM(1);
  Mt(1,2)=MM(4);
  Mt(1,3)=MM(5);
  Mt(2,1)=MM(4);
  Mt(2,2)=MM(2);
  Mt(2,3)=MM(6);
  Mt(3,1)=MM(5);
  Mt(3,2)=MM(6);
  Mt(3,3)=MM(3);
  g=9.8;
  
  zn=1 + mt' *inv(Mt) *mt;
  cc=(g^2 * mt' *Mt *mt)/zn;
  M =Mt.*(g^2/zn);
  m=mt.*(g^2/zn);
  
  Sf=chol(M);
  bf=1/2 *inv(Sf)'*m;
  
  fz=Sf*c'+bf;
  
  figure;
  plot(fz(1,:),'r');
  hold on;
  g=0;
  plot(fz(2,:),'g');
  plot(fz(3,:),'b');
end


