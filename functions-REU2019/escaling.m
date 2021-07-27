function escaling(a,b)
close all
x=-1:0.1:1;
y=abs(((b-a)*x)./((b-a)*x+(a+b)));
plot(x,y)

figure
z=a:0.1:b;
w=abs(2*z./(2*z-(a+b)));
plot(z,w)