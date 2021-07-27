function fun3(k)
close all
epsi=0:1:k-1;
y1=epsi+1+16*(epsi+1).^2;
y2=(64/3)*(k-epsi).^3.*(epsi+1).^2;

plot(epsi,y1, epsi, y2)
