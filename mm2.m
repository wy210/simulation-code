% clear;
% a11=unifrnd (0.2,0.5);
% a12=unifrnd (0.2,0.5);
% a15=1.09-a11-a12;
% 
% a21=unifrnd (0.2,0.5);
% a22=unifrnd (0.2,0.5);
% a23=1.09-a21-a22;
% 
% a32=unifrnd (0.2,0.5);
% a33=unifrnd (0.2,0.5);
% a34=1.09-a32-a33;
% 
% a43=unifrnd (0.2,0.5);
% a44=unifrnd (0.2,0.5);
% a48=1.09-a43-a44;
% 
% a51=unifrnd (0.2,0.5);
% a55=unifrnd (0.2,0.5);
% a56=1.09-a51-a55;
% 
% a65=unifrnd (0.2,0.5);
% a66=unifrnd (0.2,0.5);
% a67=1.09-a65-a66;
% 
% a76=unifrnd (0.2,0.5);
% a77=unifrnd (0.2,0.5);
% a78=1.09-a76-a77;
% 
% a87=unifrnd (0.2,0.5);
% a88=unifrnd (0.2,0.5);
% a89=1.09-a87-a88;
% 
% a98=unifrnd (0.2,0.5);
% a99=unifrnd (0.2,0.5);
% a910=1.09-a98-a99;
% 
% a109=unifrnd (0.2,1);
% a1010=1.09-a109;
% 
% A=[a11 a12 0 0 a15 0 0 0 0 0;...,
%    a21 a22 a23 0 0 0 0 0 0 0;...,
%    0 a32 a33 a34 0 0 0 0 0 0;...,
%    0 0 a43 a44 0 0 0 a48 0 0;...,
%    a51 0 0 0 a55 a56 0 0 0 0;...,
%    0 0 0 0 a65 a66 a67 0 0 0;...,
%    0 0 0 0 0 a76 a77 a78 0 0;...,
%    0 0 0 0 0 0 a87 a88 a89 0;...,
%    0 0 0 0 0 0 0 a98 a99 a910;...,
%    0 0 0 0 0 0 0 0 a109 a1010];  %system matrix
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

load A
TT=100;
Pt=cell(1,TT); Kt=cell(1,TT);
Pt{1}=10*eye(10);

for t=1:TT-1
    K=ones(10).*SP;
    Pi=C*A*Pt{t}*A'*C'+C*W*C'+V;   % Intermediate variable
    E=C*A*Pt{t}*A'+C*W;   % Intermediate variable
    lt0=100000;
    lt=9999;
     while lt0>lt
         lt0=lt;
     Gradient=((Pi*(SP.*K)'-E).*SP')';
     K=K-0.001*Gradient;
     P=(A-K*C*A)*Pt{t}*(A-K*C*A)'+(eye(10)-K*C)*W*(eye(10)-K*C)'+K*V*K';
      lt=trace(P);
     end
    Kt{t}=K;
    Pt{t+1}=P;
end
        xt=cell(1,TT); hxt=cell(1,TT); yt=cell(1,TT);
        xt{1}=normrnd(0,10,10,1); hxt{1}=zeros(10,1); yt{1}=C*xt{1}+normrnd(0,1.5,10,1);
 for t=1:TT-2
     xt{t+1}=A*xt{t}+normrnd(0,1.5,10,1);
     yt{t+1}=C*xt{t+1}+normrnd(0,1.5,10,1);
     hxt{t+1}=A*hxt{t}+Kt{t+1}*(yt{t+1}-C*A*hxt{t});
 end
  
 
 
figure(1);
 xi=0:74;
 for j=1:10
  yy=[];
  for i=1:75
          yy(i)=xt{i}(j);
  end
   plot(xi,yy,'LineWidth',1);
          hold on;
 end
   
figure(2);
 xi=0:74;
 for j=1:10
  yy=[];
  for i=1:75
          yy(i)=hxt{i}(j);
  end
   plot(xi,yy,'LineWidth',1);
          hold on;
 end
 
 figure(3);
 xi=0:74;
 for j=1:10
  yy=[];
  for i=1:75
          yy(i)=xt{i}(j)-hxt{i}(j);
  end
   plot(xi,yy,'LineWidth',1);
          hold on;
 end
  
figure(4);
 xi=0:74;
  yy=[];
  for i=1:75
          yy(i)=trace(Pt{i});
  end
   plot(xi,yy,'LineWidth',1);
          hold on;

   
%ylabel('$t$','Interpreter','LaTex');
%legend('$e_{t}^{1}$','$e_{t}^{2}$','$e_{t}^{3}$','$e_{t}^{4}$','$e_{t}^{5}$','$e_{t}^{6}$','$e_{t}^{7}$','$e_{t}^{8}$','$e_{t}^{9}$','$e_{t}^{10}$','Interpreter','LaTex');

