 clc
 clear all
 close all
 
 %% Generisanje odbiraka 
 
 % K1 
 M1=[-4;4];
 S1=[2 -0.5; -0.5 2];
 
 % K2
 M2=[4;8];
 S2=[0.9 0.7; 0.7 0.9];
 
 % d2 krive
 d2=[1 4 9];
 
 % tacke
 Nx = 40; %33;
 Ny = 40; %33;
 N = Nx*Ny; % 500
 x_min = -7; x_max = 8;
 y_min = 0; y_max = 11;
 x = linspace(x_min,x_max,Nx);  
 y = linspace(y_min,y_max,Ny); 
 delta_x = x(2)-x(1);
 delta_y = y(2)-y(1); 

X1 = mvnrnd(M1,S1,N);
X2 = mvnrnd(M2,S2,N);

% X1 = X1'; 
% X2 = X2';

figure(1)
plot(X1(:,1),X1(:,2),'x'); hold on;
plot(X2(:,1),X2(:,2),'x'); 
title("Generisani odbrici");
xlabel("x1"); ylabel("x2");
legend("prva klasa", "druga klasa", 'Location','NorthWest')

% GMT
fx1 = zeros(size(X1,1),1);
fx2 = zeros(size(X2,1),1);
for i=1:size(X1,1)
    fx1(i) = 1/(2*pi*sqrt(det(S1)))*exp(-1/2*(X1(i,:)-M1')*S1^(-1)*(X1(i,:)-M1')');
    fx2(i) = 1/(2*pi*sqrt(det(S2)))*exp(-1/2*(X2(i,:)-M2')*S2^(-1)*(X2(i,:)-M2')');
end
%fx1 = 1/(2*pi*sqrt(det(S1)))*exp(-1/2*(X1-repmat(M1',[size(X1,1),1]))*S1^(-1)*(X1-repmat(M1',[size(X1,1),1]))');
fx1_1 = 1./(2*pi*sqrt(det(S1)))*exp(-1/2);
fx1_2 = 1./(2*pi*sqrt(det(S1)))*exp(-2/2);
fx1_3 = 1./(2*pi*sqrt(det(S1)))*exp(-3/2);
fx1_4 = 1./(2*pi*sqrt(det(S1)))*exp(-4/2);

%fx2 = 1./(2*pi*sqrt(det(S2)))*exp(-1/2*(X2'-M2)'.*S2^(-1).*(X2'-M2));
fx2_1 = 1./(2*pi*sqrt(det(S2)))*exp(-1/2);
fx2_2 = 1./(2*pi*sqrt(det(S2)))*exp(-2/2);
fx2_3 = 1./(2*pi*sqrt(det(S2)))*exp(-3/2);
fx2_4 = 1./(2*pi*sqrt(det(S2)))*exp(-4/2);

E = 0.05;
X1_1 = X1(abs(fx1-fx1_1)<E*fx1,:);
X1_2 = X1(abs(fx1-fx1_2)<E*fx1,:);
X1_3 = X1(abs(fx1-fx1_3)<E*fx1,:);
X1_4 = X1(abs(fx1-fx1_4)<E*fx1,:);

X2_1 = X2(abs(fx2-fx2_1)<E*fx2,:);
X2_2 = X2(abs(fx2-fx2_2)<E*fx2,:);
X2_3 = X2(abs(fx2-fx2_3)<E*fx2,:);
X2_4 = X2(abs(fx2-fx2_4)<E*fx2,:);

figure(2)
plot3(X1(:,1),X1(:,2),fx1,'b.'); hold on;
plot3(X2(:,1),X2(:,2),fx2,'b.'); hold on;
plot3(X1_1(:,1),X1_1(:,2),fx1_1*ones(size(X1_1,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X1_2(:,1),X1_2(:,2),fx1_2*ones(size(X1_2,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X1_3(:,1),X1_3(:,2),fx1_3*ones(size(X1_3,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X1_4(:,1),X1_4(:,2),fx1_4*ones(size(X1_4,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X2_1(:,1),X2_1(:,2),fx2_1*ones(size(X2_1,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X2_2(:,1),X2_2(:,2),fx2_2*ones(size(X2_2,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X2_3(:,1),X2_3(:,2),fx2_3*ones(size(X2_3,1)),'g*'); hold on;  %,'LineWidth',2); hold on;
plot3(X2_4(:,1),X2_4(:,2),fx2_4*ones(size(X2_4,1)),'g*'); hold off;  %,'LineWidth',2); hold on;
title('GMT')
grid on

% plot(fx1_2(:,1),fx1_2(:,2),'b','LineWidth',2); hold on;
% plot(fx1_3(:,1),fx1_3(:,2),'b','LineWidth',2); hold on;
% plot(fx1_4(:,1),fx1_4(:,2),'b','LineWidth',2); hold on;
% plot(fx2_1(:,1),fx2_1(:,2),'b','LineWidth',2); hold on;
% plot(fx2_2(:,1),fx2_2(:,2),'b','LineWidth',2); hold on;
% plot(fx2_3(:,1),fx2_3(:,2),'b','LineWidth',2); hold on;
% plot(fx2_4(:,1),fx2_4(:,2),'b','LineWidth',2); hold on;
 
 %% BAJESOV TEST MINIMALNE GRESKE
 
 % inicijalizacija potrebnih promenljivih
 f1 = zeros(Nx,Ny);
 f2 = zeros(Nx,Ny);
 h_bajes = zeros(Nx,Ny);
 X1_bajes = [];
 X2_bajes = [];
 Eps1_bajes = 0;
 Eps2_bajes = 0;
 
 for i=1:length(x)
     for j=1:length(y)
         X=[x(i);y(j)];
         f1(i,j)=1/(2*pi*det(S1)^0.5)*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
         f2(i,j)=1/(2*pi*det(S2)^0.5)*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));
         h_bajes(i,j)=-log(f1(i,j))+ log(f2(i,j)); 
         if h_bajes(i,j)<0
             X1_bajes = [X1_bajes, X];
             Eps2_bajes=Eps2_bajes+f2(i,j)*delta_x*delta_y;
         else
             X2_bajes = [X2_bajes, X];
             Eps1_bajes=Eps1_bajes+f1(i,j)*delta_x*delta_y;
         end
     end
 end
 
 figure(5);
 mesh(x,y,f1'); hold on; mesh(x,y,f2'); hold on;
 xlabel("x"); ylabel("y");
 contour(x,y,h_bajes',[0 1 -1],'LineWidth',3); hold off;
 
 Eps1_bajes=0;
 Eps2_bajes=0;
 for i =1:length(x)
     for j=1:length(y)
         if h_bajes(i,j)<0 
             Eps2_bajes=Eps2_bajes+f2(i,j)*delta_x*delta_y;
         else
             Eps1_bajes=Eps1_bajes+f1(i,j)*delta_x*delta_y;
             
         end
     end
 end
 
 % d2 krive
 prag1=max(max(f1))*exp(-0.5*d2);
 prag2=max(max(f2))*exp(-0.5*d2);

 figure(3)
 plot(X1_bajes(1,:),X1_bajes(2,:),'r.'); hold on; 
 plot(X2_bajes(1,:),X2_bajes(2,:),'b.'); hold on; 
 contour(x,y,h_bajes',[0 1 -1]); hold on;
 contour(x,y,f1',prag1,'Linewidth',3); hold on;
 contour(x,y,f2',prag2,'Linewidth',3); hold on;
 plot(X1(:,1),X1(:,2),'rx'); hold on;
 plot(X2(:,1),X2(:,2),'bx'); 
 xlabel("x"); ylabel("y");
 legend('K1','K2','klasifikaciona linija','GMT','Location','SouthEast')
% legend('K1','d2 krive K1','K2','d2 krive K2','klasifikaciona linija',...
%     'Location','SouthEast')
 title('Klasifikacija Basejovim testom minimalne verovatnoce greske')
 xlim([-7,8]); ylim([0,11])

 %% KLASIFIKATOR DISTANCE
 
 h_dist = zeros(Nx,Ny);
 X1_dist = [];
 X2_dist = [];
 Eps1_dist = 0;
 Eps2_dist = 0;
 
 
 for i=1:length(x)
     for j=1:length(y)
         X=[x(i);y(j)];
         h_dist(i,j)= (X-M2)'*(X-M2)-(X-M1)'*(X-M1);
         if h_dist(i,j)>=0
             X1_dist = [X1_dist, X];
             Eps2_dist=Eps2_dist+f2(i,j)*delta_x*delta_y;
         else
             X2_dist = [X2_dist, X];
             Eps1_dist=Eps1_dist+f1(i,j)*delta_x*delta_y;
         end
     end
 end
 
 figure(4)
 plot(X1_dist(1,:),X1_dist(2,:),'r.'); hold on; 
 plot(X2_dist(1,:),X2_dist(2,:),'b.'); hold on; 
 contour(x,y,h_dist',[0 1 -1]); hold on;
 contour(x,y,f1',prag1,'Linewidth',3); hold on;
 contour(x,y,f2',prag2,'Linewidth',3); hold on;
 plot(X1(:,1),X1(:,2),'rx'); hold on;
 plot(X2(:,1),X2(:,2),'bx');
 xlabel("x"); ylabel("y");
 legend('K1','K2','klasifikaciona linija','GMT','Location','SouthEast')
%  legend('K1','d2 krive K1','K2','d2 krive K2','klasifikaciona linija',...
%      'Location','SouthEast')
 title('Klasifikator distance')
 xlim([-7,8]); ylim([0,11])

%  Eps1_dist=0;
%  Eps2_dist=0;
%  for i =1:length(x)
%      for j=1:length(y)
%          if h_dist(i,j)<0 
%              Eps2_dist=Eps2_dist+f2(i,j)*delta_x*delta_y;
%          else
%              Eps1_dist=Eps1_dist+f1(i,j)*delta_x*delta_y;
%              
%          end
%      end
%  end
 
%% Wald-ov sekvencijalni test
%  
%  e_min = 1e-5;
%  e1_vector = linspace(e_min,1-e_min,N);
%  e2_vector = linspace(e_min,1-e_min,N); 
%  m1_eps1_vector = zeros(1,length(e1_vector));
%  m2_eps1_vector = zeros(1,length(e1_vector));
%  m1_eps2_vector = zeros(1,length(e2_vector));
%  m2_eps2_vector = zeros(1,length(e2_vector));
%  lm1_vector = 100*ones(20,N);
%  lm2_vector = 100*ones(20,N);
%  e_fix = 1e-5;
%  X1_wald = [];
%  X2_wald = [];
%  X_curr = [];
%  m_vector = [];
%  N1_bajes = length(X1_bajes(1,:));
%  N2_bajes = length(X2_bajes(1,:));
%  
%  for k = 1:length(e1_vector)
%      % granice za Wald-ov test
%      e1 = e1_vector(k);
%      A = (1-e1)/e_fix;
%      B = e1/(1-e_fix);
%      % mesanje vektora
%      ind1 = randperm(N1_bajes);
%      X1 = X1_bajes(:,ind1);
%      ind2 = randperm(N2_bajes);
%      X2 = X2_bajes(:,ind2);
%      % inicijalizacija promenljivih
%      lm = 1;
%      X_curr = [];
%      m1 = [];
%      m2 = [];
%      lm_curr = 100*ones(20,1);
%      for i = 1:N1_bajes
%          X=X1(:,i);
%          X_curr = [X_curr X];
%          f1 = 1/(2*pi*det(S1)^0.5)*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
%          f2 = 1/(2*pi*det(S2)^0.5)*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));
%          lm = lm*f1/f2;
%          lm_curr =[];
%          if (lm>=A)
%              X1_wald = [X1_wald X_curr];
%              m1 = [m1 i];
%              lm = 1;
%              X_curr = [];
%          end
%      end
%      
%      for j = 1:N2_bajes
%          X=X2(:,j);
%          X_curr = [X_curr X];
%          f1 = 1/(2*pi*det(S1)^0.5)*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
%          f2 = 1/(2*pi*det(S2)^0.5)*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));
%          lm = lm*f1/f2;
%          lm_vector(i,k) = lm;
%          if (lm<=B) % ne bi trebalo nikad da se desi 
%              X2_wald = [X2_wald X_curr];
%              m2 = [m2 j];
%              lm = 1;
%              X_curr = [];
%          end
%      end
%      
%      
%      m1_eps1_vector(k) = sum(m1)/length(m1);
%      
%      m2_eps1_vector(k) = sum(m2)/length(m2);
%      
%  end
%  
%  
%   for k = 1:length(e2_vector)
%      % granice za Wald-ov test
%      e2 = e2_vector(k);
%      A = (1-e_fix)/e2;
%      B = e_fix/(1-e2);
%      % mesanje vektora
%      ind1 = randperm(N1_bajes);
%      X1 = X1_bajes(:,ind1);
%      ind2 = randperm(N2_bajes);
%      X2 = X2_bajes(:,ind2);
%  
%      % inicijalizacija promenljivih
%      lm = 1;
%      X_curr = [];
%      m1 = [];
%      m2 = [];
%      for i = 1:N1_bajes
%          X=X1(:,i);
%          X_curr = [X_curr X];
%          f1 = 1/(2*pi*det(S1)^0.5)*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
%          f2 = 1/(2*pi*det(S2)^0.5)*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));
%          lm = lm*f1/f2;
%          if (lm>=A)
%              X1_wald = [X1_wald X_curr];
%              m1 = [m1 i];
%              lm = 1;
%              X_curr = [];
%          end
%      end
%      
%      for j = 1:N2_bajes
%          X=X2(:,j);
%          X_curr = [X_curr X];
%          f1 = 1/(2*pi*det(S1)^0.5)*exp(-0.5*(X-M1)'*inv(S1)*(X-M1));
%          f2 = 1/(2*pi*det(S2)^0.5)*exp(-0.5*(X-M2)'*inv(S2)*(X-M2));
%          lm = lm*f1/f2;
%          if (lm<=B) % ne bi trebalo nikad da se desi 
%              X2_wald = [X2_wald X_curr];
%              m2 = [m2 j];
%              lm = 1;
%              X_curr = [];
%          end
%      end
%      
%      m1_eps2_vector(k) = sum(m1)/length(m1);
%      
%      m2_eps2_vector(k) = sum(m2)/length(m2);
%      
%  end
%  
%  
%  %%
%  figure(5)
%  plot(e1_vector,log(m1_eps1_vector))
%  title('Log zavisnost m1 u zavisnosti od greske prvog tipa')
%  figure(5)
%  plot(e1_vector,log(m2_eps1_vector))
%  title('Log zavisnost m2 u zavisnosti od greske prvog tipa')
%  figure(6)
%  plot(e2_vector,log(m1_eps2_vector))
%  title('Log zavisnost m1 u zavisnosti od greske drugog tipa')
%  figure(7)
%  plot(e2_vector,log(m2_eps2_vector))
%  title('Log zavisnost m2 u zavisnosti od greske drugog tipa')
%  
