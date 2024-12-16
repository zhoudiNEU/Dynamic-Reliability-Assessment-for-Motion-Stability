clc
clear
ALLA=[];
for iii=-0.1:0.001:0.05
    iii
A=readmatrix('C:\Users\86159\Desktop\Data\优化\90\90+80+40+10.csv');
B=A(A>5);
mm = mean(B+iii);
% biaozhuncha1
ee=cov(B);
% sanjiju
ss= mean( (B-mm).^3 );
% sijiju
sss= mean( (B-mm).^4 );


muX=[mm;9];
cvX=[ee;0.1];
sigmaX = cvX.*muX;
csX=[ss;0.2];
ckX=[sss;0.2];
m=4;
gX=[1;-1];
gXX=zeros(2);
muG=muX(1)-muX(2);
m2=sigmaX.^2;m3=csX.*sigmaX.^3;m4=ckX.*sigmaX.^4;
t42=m4-3*m2.^2;dg=diag(gXX);gX2=gX.^2;gt2=gX.*m2;
muZ=muG+dg.'*m2/2;
muZ2=gX2.'*m2+gX.'.*dg.'*m3+dg.'.^2*t42/4+m2.'*gXX.^2*m2/2;
sigmaZ= sqrt(muZ2);
muZ3=gX.'.^3*m3+1.5*gX2.'.*dg.'*t42+3*gt2.'*gXX*gt2;
muZ4= gX.'.^4*t42+3*gt2.'*gX*gX.'*gt2;
csZ=muZ3/sigmaZ^3;
ckZ= muZ4/sigmaZ^4;
nuY=[1;0;1;csZ;ckZ];
t=10*ckZ-12*csZ^2-18;c0=(3*csZ^2-4*ckZ)/t;
c1=-csZ*(3+ckZ)/t;c2=(6+3*csZ^2-2*ckZ)/t;
for k=4:m-1
    nuY(k+2)=-k*(c0*nuY(k)+c1*nuY(k+1))/(1+(k+2)*c2);
end
a = fsolve(@(a)fint(a,nuY,m,10),zeros(m+1,1));
pF = integral(@(z)exp(-polyval(a,z)),-Inf,-muZ/sigmaZ);

ALLA=[ALLA,pF];
end
AALL=[1,ALLA(1:1:30)];
plot(0:10:300,AALL,'b')
