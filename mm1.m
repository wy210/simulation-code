clear;
a11=unifrnd (0.2,0.5);
a12=unifrnd (0.2,0.5);
a15=1.09-a11-a12;

a21=unifrnd (0.2,0.5);
a22=unifrnd (0.2,0.5);
a23=1.09-a21-a22;

a32=unifrnd (0.2,0.5);
a33=unifrnd (0.2,0.5);
a34=1.09-a32-a33;

a43=unifrnd (0.2,0.5);
a44=unifrnd (0.2,0.5);
a48=1.09-a43-a44;

a51=unifrnd (0.2,0.5);
a55=unifrnd (0.2,0.5);
a56=1.09-a51-a55;

a65=unifrnd (0.2,0.5);
a66=unifrnd (0.2,0.5);
a67=1.09-a65-a66;

a76=unifrnd (0.2,0.5);
a77=unifrnd (0.2,0.5);
a78=1.09-a76-a77;

a87=unifrnd (0.2,0.5);
a88=unifrnd (0.2,0.5);
a89=1.09-a87-a88;

a98=unifrnd (0.2,0.5);
a99=unifrnd (0.2,0.5);
a910=1.09-a98-a99;

a109=unifrnd (0.2,1);
a1010=1.09-a109;

A=[a11 a12 0 0 a15 0 0 0 0 0;...,
   a21 a22 a23 0 0 0 0 0 0 0;...,
   0 a32 a33 a34 0 0 0 0 0 0;...,
   0 0 a43 a44 0 0 0 a48 0 0;...,
   a51 0 0 0 a55 a56 0 0 0 0;...,
   0 0 0 0 a65 a66 a67 0 0 0;...,
   0 0 0 0 0 a76 a77 a78 0 0;...,
   0 0 0 0 0 0 a87 a88 a89 0;...,
   0 0 0 0 0 0 0 a98 a99 a910;...,
   0 0 0 0 0 0 0 0 a109 a1010];  %system matrix
% v=[1 1 1 1 1 1 1 1 1 1]*2.48;
% C=diag(v);   % Measurement matrix
W=1.5*eye(10);   % covariance of process noise
V=1.5*eye(10);   % covariance of measurement noise
SP=[1 1 0 0 1 0 0 0 0 0;1 1 1 0 0 0 0 0 0 0;...,
    0 1 1 1 0 0 0 0 0 0;0 0 1 1 0 0 0 1 0 0;...,
    1 0 0 0 1 1 0 0 0 0;0 0 0 0 1 1 1 0 0 0;...,
    0 0 0 0 0 1 1 1 0 0;0 0 0 0 0 0 1 1 1 0;...,
    0 0 0 0 0 0 0 1 1 1;0 0 0 0 0 0 0 0 1 1];
C=2.6*SP;

save A;

da=0.0001;
EI=eye(10);

FF=@ (x,y) (A-(SP.*x)*C*A)*y*(A-(SP.*x)*C*A)';
HH=@ (x,e) FF(x,e*e'+da*EI)-e*e'+da*EI;
Gr=@ (x,y) (-C*A*y*(A-(x.*SP)*C*A)')'.*SP;
Ge=@ (x,e) ((A-(SP.*x)*C*A)*(A-(SP.*x)*C*A)'-EI)*e;

for i=1:200
ee=normrnd(0.5,6,10,10);
pp=10000000*eye(10);
pp1=10000001*eye(10);
while  max(eig(pp))<max(eig(pp1))
       pp1=pp;
yy=ee*ee'++da*EI;
  xx=0.1*ones(10).*SP;
    lt0=100000;
    lt=9999;
     while lt0>lt && lt0-lt>0.0005
         lt0=lt;
     gg=Gr(xx,yy);
     gg=gg/norm(gg);
    xx=xx-0.015*gg;
     hh=HH(xx,ee);
      lt=trace(hh);
     end
   pp=hh;
if max(eig(hh))<0
    disp('true');
    return;
end
ge=Ge(xx,ee);
ge=ge/norm(ge);
ee=ee-0.001*ge;
end
end


