clc 
clear all
close all

%% Generisanje odbiraka klasa

N = 1000;

M1 = [0;3.5];
s1_1 = 1;
s2_1 = 0.8;
ro_1 = 0.6;
S1 = [s1_1^2 ro_1*s1_1*s2_1; ro_1*s1_1*s2_1 s2_1^2];
X1 = mvnrnd(M1,S1,N);

M2 = [-2;-1];
s1_2 = 1;
s2_2 = 1;
ro_2 = -0.4;
S2 = [s1_2^2 ro_2*s1_2*s2_2; ro_2*s1_2*s2_2 s2_2^2];
X2 = mvnrnd(M2,S2,N); 

X1 = X1';
X2 = X2';

figure, hold all
plot(X1(1,:),X1(2,:),'x'); 
plot(X2(1,:),X2(2,:),'x'); 
legend('klasa 1', 'klasa 2')
title('Raspored odbiraka')
xlabel('x1'); ylabel('x2')

% %% Podela na obucavajuci i test skup
% 
% ind = randperm(N);
% Y1trening = Y1(ind(1:0.7*N),:);
% Y1test = Y1(ind(0.7*N+1:N),:);
% 
% ind = randperm(N);
% Y2trening = Y2(ind(1:0.7*N),:);
% Y2test = Y2(ind(0.7*N+1:N),:);

%% projektovanje klasifikatora

V = zeros(1,2);
% sigma1 = 0;
% sigma2 = 0;

Ne_best = N+1; %bice sigurno manje od ove vrednosti 
Nes = zeros(1,100);
v = zeros(1,100);
i = 0;

for s = 0:0.01:1
    i = i+1;
    
    V = (s*S1+(1-s)*S2)^(-1)*(M2-M1);
%     sigma1 = V'*S1*V;
%     sigma2 = V'*S2*V;
    Y1 = V'*X1;
    Y2 = V'*X2;
    
%     v0min = -max(max(sqrt(Y1(1,:).^2+Y1(2,:).^2)), max(sqrt(Y2(1,:).^2+Y2(2,:).^2)));
%     v0max(1) = -min(min(sqrt(Y1(1,:).^2+Y1(2,:).^2)), min(sqrt(Y2(1,:).^2+Y2(2,:).^2)));
    v0min = -100;
    v0max = 100;
    
    Nebest = N+1;
    for v0 = v0min:0.01:v0max
        
        Ne = length(Y1(Y1>-v0));
        Ne = Ne + length(Y2(Y2<-v0));
        if (Ne<Ne_best)
            Ne_best = Ne;
            v(i) = v0;
        end
    end
    Nes(i) = Ne_best;
end

[~, indMin] = min(Nes);
s = 0+0.01*indMin;
V = (s*S1+(1-s)*S2)^(-1)*(M2-M1);
vo = v(indMin);

%%
figure(2)
plot(linspace(0,1,101), Nes)
title('Zavisnost procenjene greske od vrednosti parametra s')

V
vo

%% 

xp1 = [-5 5]; 

figure(3), hold all
plot(X1(1,:),X1(2,:),'x'); 
plot(X2(1,:),X2(2,:),'x'); 
plot(xp1,(-vo-V(1)*xp1)/V(2), 'LineWidth',2);
legend('klasa 1', 'klasa 2')
title('Raspored odbiraka')
xlabel('x1'); ylabel('x2')



%% b) Metoda zeljenog izlaza 

N1 = length(X1(1,:));
N2 = length(X2(1,:));
U = [-ones(1,N1), ones(1,N2); -X1, X2];
Gamma = [ones(N1,1); ones(N2,1)];
W = inv(U*U')*U*Gamma;
%%
figure(4), hold all
plot(X1(1,:),X1(2,:),'x'); 
plot(X2(1,:),X2(2,:),'x'); 
plot(xp1,(-W(1)-W(2)*xp1)/W(3), 'LineWidth',2);
legend('klasa 1', 'klasa 2')
title('Klasifikacija metodom zeljenog izlaza')
xlabel('x1'); ylabel('x2')

Gamma = [2*ones(N1,1); ones(N2,1)];
W = inv(U*U')*U*Gamma;
%%
figure(5), hold all
plot(X1(1,:),X1(2,:),'x'); 
plot(X2(1,:),X2(2,:),'x'); 
plot(xp1,(-W(1)-W(2)*xp1)/W(3), 'LineWidth',2);
legend('klasa 1', 'klasa 2')
title('Klasifikacija kad je znacajnija prva klasa')
xlabel('x1'); ylabel('x2')

Gamma = [ones(N1,1); 2*ones(N2,1)];
W = inv(U*U')*U*Gamma;

figure(6), hold all
plot(X1(1,:),X1(2,:),'x'); 
plot(X2(1,:),X2(2,:),'x'); 
plot(xp1,(-W(1)-W(2)*xp1)/W(3), 'LineWidth',2);
legend('klasa 1', 'klasa 2')
title('Klasifikacija kad je znacajnija druga klasa')
xlabel('x1'); ylabel('x2')

%% 2. Kvadratni klasifikator

N = 1000; 

% Klasa 1 (polukrug)
Mx = [8;8]; 
Rx = 3*rand(1,N);
alfax = pi*rand(1,N); 
X = [Rx.*cos(alfax); Rx.*sin(alfax)] + Mx*ones(1,N); 

% Klasa 2 (poluprsten)
My = [8;8]; 
R_u = 5; 
R_d = 2; 
Ry = R_u + R_d*rand(1,N);
alfay = pi*rand(1,N); 
Y = [Ry.*cos(alfay); Ry.*sin(alfay)] + My*ones(1,N); 
figure(7);
plot(X(1,:),X(2,:),'bx'); hold on
plot(Y(1,:),Y(2,:),'rx');
title('Raspored odbiraka')

Gamma = [ones(2*N,1)]; 
Z = [-1*ones(1,N), 1*ones(1,N); -X, Y; -X(1,:).^2, Y(1,:).^2; ...
    2*(-X(1,:)).*(X(2,:)), 2*(Y(1,:)).*(Y(2,:)); -X(2,:).^2, Y(2,:).^2];
W = (Z*Z')^(-1)*Z*Gamma;

v0 = W(1);
v = W(2:3);
Q = [W(4) W(5); W(5) W(6)];

% iscrtavanje
syms xp yp 
[xp,yp,~,~] = solve(v0+xp*v(1)+yp*v(2)+xp^2*Q(1)+2*xp*yp*Q(2)+yp^2*Q(4),xp,yp,'returnConditions',true);

z = -10:0.001:25;
xp = eval(xp);
xp = [xp(1,:), fliplr(xp(2,:))];
xp1 = xp(imag(xp)==0); 
yp = eval(yp);
yp = [yp(1,:), fliplr(yp(2,:))];
yp1 = yp(imag(xp)==0);

figure(8)
plot(Y(1,:),Y(2,:),'rx'); hold on
plot(X(1,:),X(2,:),'bx');
plot(xp1,yp1,'k','LineWidth',1);
ylim([8 15])
title('Klasifikacija kvadratnim klasifikatorom')














