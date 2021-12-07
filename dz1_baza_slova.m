clear all
close all
clc
%%
% Definisemo broj klasa je 5
Nk = 5;
% Definisemo broj obelezja
Nobelezja = 5;
% Ucitavamo slova A
x=dir('bazaA*.bmp');
P1 = zeros(Nobelezja,max(size(x)));
% one koje zelimo da plotujemo oznacicemo 1 u nizu a
a = zeros(1,max(size(x)));
a([105,118]) = 1; 
for i=1:max(size(x))
    X=imread(x(i).name); % name=['bazaA' num2str(i,'%03d') '.bmp'];
    P1(:,i)=Obelezja(X,a(i));
end


% Ucitavamo slova E
x=dir('bazaE*.bmp');
P2 = zeros(Nobelezja,max(size(x)));
% one koje zelimo da plotujemo oznacicemo 1 u nizu a
a = zeros(1,max(size(x)));
a([105,109,110,113]) = 1; 
for i=1:max(size(x))
    X=imread(x(i).name); 
    P2(:,i)=Obelezja(X,a(i));
end


% Ucitavamo slova I
x=dir('bazaI*.bmp');
P3 = zeros(Nobelezja,max(size(x)));
% one koje zelimo da plotujemo oznacicemo 1 u nizu a
a = zeros(1,max(size(x)));
a([101,106,109,117]) = 1; 
for i=1:max(size(x))
    X=imread(x(i).name); 
    P3(:,i)=Obelezja(X,a(i));
end


% Ucitavamo slova O
x=dir('bazaO*.bmp');
P4 = zeros(Nobelezja,max(size(x)));
for i=1:max(size(x))
    X=imread(x(i).name); 
    P4(:,i)=Obelezja(X,0);
end


% Ucitavamo slova U
x=dir('bazaU*.bmp');
P5 = zeros(Nobelezja,max(size(x)));
for i=1:max(size(x))
    X=imread(x(i).name); 
    P5(:,i)=Obelezja(X,0);
end

%% Plotovi obelezja

figure(1);
plot(P1(1,:),P1(2,:),'ro'); hold on;
plot(P2(1,:),P2(2,:),'bx');
plot(P3(1,:),P3(2,:),'gv');
plot(P4(1,:),P4(2,:),'kx');
plot(P5(1,:),P3(2,:),'mv'); hold off;
legend('A','E','I','O','U','Location','SouthEast');
xlabel('prvo obelezje'); ylabel('drugo obelezje')

figure(2);
plot(P1(3,:),P1(4,:),'ro'); hold on;
plot(P2(3,:),P2(4,:),'bx');
plot(P3(3,:),P3(4,:),'gv');
plot(P4(3,:),P4(4,:),'kx');
plot(P5(3,:),P5(4,:),'mv'); hold off;
legend('A','E','I','O','U','Location','SouthEast');
xlabel('trece obelezje'); ylabel('cetvrto obelezje')

figure(3);
plot(P1(4,:),P1(5,:),'ro'); hold on;
plot(P2(4,:),P2(5,:),'bx');
plot(P3(4,:),P3(5,:),'gv');
plot(P4(4,:),P4(5,:),'kx');
plot(P5(4,:),P5(5,:),'mv'); hold off;
legend('A','E','I','O','U','Location','SouthEast');
xlabel('cetvrto obelezje'); ylabel('peto obelezje')

%%
N=120;
No=100; % podela na O:T skup obicno ide u odnosu 60:40/70:30/80:20
O1=P1(:,1:No);T1=P1(:,No+1:N);
O2=P2(:,1:No);T2=P2(:,No+1:N);
O3=P3(:,1:No);T3=P3(:,No+1:N);
O4=P4(:,1:No);T4=P4(:,No+1:N);
O5=P5(:,1:No);T5=P5(:,No+1:N);

%P1=P2=P3=P4=5, raspodele su normalne

M1=mean(O1,2); S1=cov(O1'); % mean(O1')';
M2=mean(O2,2); S2=cov(O2');
M3=mean(O3,2); S3=cov(O3');
M4=mean(O4,2); S4=cov(O4');
M5=mean(O5,2); S5=cov(O5');


% Testiranje klasifikatora: za svaku klasu, testiramo svaki odbirak
% i belezimo da li je odluka klasifikatora ispravna ili ne

Mk=zeros(Nk); % konfuziona matrica
for kl=1:5
    disp(num2str(kl));
    if kl==1
        T=T1;
    elseif kl==2
        T=T2;
    elseif kl==3
        T=T3;
    elseif kl==4
        T=T4;
    else
        T=T5;
    end
    for i=1:N-No %           
        p=T(:,i);
        f1=1/((2*pi)^(Nobelezja/2))/det(S1)^0.5*exp(-0.5*(p-M1)'*inv(S1)*(p-M1));
        f2=1/((2*pi)^(Nobelezja/2))/det(S2)^0.5*exp(-0.5*(p-M2)'*inv(S2)*(p-M2));
        f3=1/((2*pi)^(Nobelezja/2))/det(S3)^0.5*exp(-0.5*(p-M3)'*inv(S3)*(p-M3));
        f4=1/((2*pi)^(Nobelezja/2))/det(S4)^0.5*exp(-0.5*(p-M4)'*inv(S4)*(p-M4));
        f5=1/((2*pi)^(Nobelezja/2))/det(S5)^0.5*exp(-0.5*(p-M5)'*inv(S5)*(p-M5));
        m=max([f1 f2 f3 f4 f5]); % maksimum izmedju svih pet daje odluku
        if m==f1 % klasifikator doneo odluku da je klasa 1
            Mk(kl,1)=Mk(kl,1)+1;
            disp('Odluka je 1');
        elseif m==f2 % klasifikator doneo odluku da je klasa 2
            Mk(kl,2)=Mk(kl,2)+1;
            disp('Odluka je 2');
        elseif m==f3 % klasifikator doneo odluku da je klasa 3
            Mk(kl,3)=Mk(kl,3)+1;
            disp('Odluka je 3');
        elseif m==f4 % klasifikator doneo odluku da je klasa 4
            Mk(kl,4)=Mk(kl,4)+1;
            disp('Odluka je 4');
        else % klasifikator doneo odluku da je klasa 5
            Mk(kl,5)=Mk(kl,5)+1;
            disp('Odluka je 5');
        end
    end
end
Mk
disp(['Greska iznosi: ' num2str((sum(sum(Mk))-trace(Mk))/sum(sum(Mk)))]);
