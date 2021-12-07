clc 
clear all
close all
 
N = 500; % br odbiraka

% pp unif raspodelu u krugu 

% Klasa 1 (ispunjen polukrug)
Mx = [8;9]; % centar
% poluprecnik je 2
Rx = 1.8*rand(1,N); % vektor radijusa
alfax = 2*pi*rand(1,N); % vektor uglova
X1 = [Rx.*cos(alfax); Rx.*sin(alfax)] + Mx*ones(1,N); % Mx je 2x1 pa moramo da ponovimo N puta 
figure(1)
plot(X1(1,:),X1(2,:),'bx');
title('Raspored odbiraka')

% Klasa 2 (poluprsten)
My = [8;8]; % centar
R_u = 5; % unutrasnji poluprecnik
R_d = 3; % sirina prstena
Ry = R_u + R_d*rand(1,N); 
alfay = pi*rand(1,N); 
X2 = [Ry.*cos(alfay); Ry.*sin(alfay)] + My*ones(1,N); % Mx je 2x1 pa moramo da ponovimo N puta 
figure(1); hold on
plot(X2(1,:),X2(2,:),'ro');
legend('klasa 1', 'klasa 2')
X_odb = [X1 X2];

 %% Klasterizacija
 
 % broj klasa
L = 2;

% Ulaz
X = [X1 X2; ones(1,N), 2*ones(1,N)];

Npocetni = 1;
Ycurr_vector = ones(1,2*N);
Ycurr_vector(X(2,:)>9) = 2;

figure, 
for i = 1:L
    plot(X(1,Ycurr_vector(1,:)==i),X(2,Ycurr_vector(1,:)==i),'x'); hold on;
end
legend('klasa 1', 'klasa 2')
title('Inicijalna klasterizacija')
xlabel('x_1'); ylabel('x_2')
 


Yfinal = zeros(Npocetni,L*N);
P = zeros(1,Npocetni);
cnt = zeros(1,Npocetni)
for i = 1:Npocetni
    Ycurr = Ycurr_vector; %((i-1)*L*N+1:i*L*N);
    [Yfinal(i,:),~,P(i),cnt(i)] = normal_decomposition(X,Ycurr,L,N,10000); % ML_clustering(X,Ycurr,L,N,0);
end

P_avg = mean(P);
cnt_avg = mean(cnt);

figure, 
for i = 1:L
    plot(X(1,Yfinal(1,:)==i),X(2,Yfinal(1,:)==i),'x'); hold on;
end
legend('klasa 1', 'klasa 2')
title('Klasterizovani odbirci')
xlabel('x_1'); ylabel('x_2')
 
 
 
 