clc
clear all 
close all 

N = 500;

%% Generisanje odbiraka klasa

M1 = [2;3.5];
s1_1 = 1;
s2_1 = 2;
ro_1 = 0.6;
S1 = [s1_1^2 ro_1*s1_1*s2_1; ro_1*s1_1*s2_1 s2_1^2];
X1 = mvnrnd(M1,S1,N);

M2 = [-4;-5];
s1_2 = 1.5;
s2_2 = 1;
ro_2 = -0.4;
S2 = [s1_2^2 ro_2*s1_2*s2_2; ro_2*s1_2*s2_2 s2_2^2];
X2 = mvnrnd(M2,S2,N); 

M3 = [-4;4];
s1_3 = 0.9;
s2_3 = 1.2;
ro_3 = -0.6;
S3 = [s1_3^2 ro_3*s1_3*s2_3; ro_3*s1_3*s2_3 s2_3^2];
X3 = mvnrnd(M3,S3,N);

M4 = [2;-5];
s1_4 = 0.7;
s2_4 = 2;
ro_4 = 0.7;
S4 = [s1_4^2 ro_4*s1_4*s2_4; ro_4*s1_4*s2_4 s2_4^2];
X4 = mvnrnd(M4,S4,N); 

X1 = X1';
X2 = X2';
X3 = X3';
X4 = X4';
%%
figure, hold all
plot(X1(1,:),X1(2,:),'x'); 
plot(X2(1,:),X2(2,:),'x'); 
plot(X3(1,:),X3(2,:),'x'); 
plot(X4(1,:),X4(2,:),'x'); 
legend('klasa 1', 'klasa 2', 'klasa 3', 'klasa 4','Location','SouthWest')
title('Raspored odbiraka')
xlabel('x1'); ylabel('x2')


%% Ispitivanje zavisnosti osetljivosti od pocetne klasterizacije

% broj klasa
L = 4;

% Ulaz
X = [X1 X2 X3 X4; ones(1,N), 2*ones(1,N) 3*ones(1,N) 4*ones(1,N)];

Npocetni = 20;
Ycurr_vector = randi([1,L],1,L*N*Npocetni);
Yfinal = zeros(Npocetni,L*N);
P = zeros(1,Npocetni);
cntv = zeros(1,Npocetni);
for i = 1:Npocetni
    Ycurr = Ycurr_vector((i-1)*L*N+1:i*L*N);
    [Yfinal(i,:),j,P(i),cntv(i)] = cmean(X,Ycurr,L,N,0);
end

figure
plot(linspace(1,20,20),P,'.')
title('Zavisnost osetljivosti od pocetne klasterizacije')
%% iscrtavanje jednog rezultata klasterizacije 
% figure, 
% 
% for i = 1:L
%     plot(X(1,Yfinal(9,:)==i),X(2,Yfinal(9,:)==i),'x'); hold on;
% end
% 
% legend('klasa 1', 'klasa 2', 'klasa 3', 'klasa 4')
% title('Klasterizovani odbirci')
% xlabel('x1'); ylabel('x2')

%% Ispitivanje zavisnosti osetljivosti od maksimalnog broja iteracija

% broj klasa
L = 4;

% ulaz
X = [X1 X2 X3 X4; ones(1,N), 2*ones(1,N) 3*ones(1,N) 4*ones(1,N)];

% pocetna klasterizacija
Ycurr = randi([1,L],1,L*N);
figure
title('Inicijalna klasterizacija')
for i = 1:L
    plot(X(1,Ycurr(:)==i),X(2,Ycurr(:)==i),'x'); hold on;
end
  
Nit_vector = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 50 80 100];
Nit = length(Nit_vector);
Yfinal = zeros(Nit,L*N);
Pit = zeros(1,Nit);
for i = 1:Nit
    [Yfinal(i,:),j,Pit(i),~] = cmean(X,Ycurr,L,N,Nit_vector(i));
end


figure
plot(Nit_vector,Pit)
title('Zavisnost osetljivosti od maksimalnog broja iteracija')


%% Ispitivanje zavisnosti osetljivosti od broja klasa

NL = 10;
Nit = 10;
cntv2 = zeros(Nit,NL);
for L = 1:NL
    Ycurr = randi([1,L],1,4*N);
    for i = 1:Nit
        [Yfinal,~,~,cntv2(i,L)] = cmean(X,Ycurr,L,N,0);
    end
    figure
    for i = 1:L
        plot(X(1,Yfinal(:)==i),X(2,Yfinal(:)==i),'x'); hold on;
    end
    hold off;
    xlabel('x_1'); ylabel('x_2')
    
end

cntv_avg = mean(cntv2);

figure
stem(cntv_avg)
title('Zavisnost prosecnog broja iteracija od broja klasa')
xlabel('L'); ylabel('broj iteracija')

