function P = Obelezja(X,a)


Y=double(X);
% figure(1);imshow(X)
% figure(2);imshow(Y/255)

% binarizacija
T=0.8;
Y(Y<T*max(max(Y)))=0;
Y(Y>=T*max(max(Y)))=255;

% predobrada slike
[nr,nc]=size(Y);
% otklanjanje okvira
poc=1;
while (poc<nr) && (sum(Y(poc,:))/nc<200)
    poc=poc+1;
end

kraj=nr;
while (kraj>1) && (sum(Y(kraj,:))/nc<200)
    kraj=kraj-1;
end

levo=1;
while (levo<nc) && (sum(Y(:,levo))/nr<200)
    levo=levo+1;
end

desno=nc;
while (desno>1) && (sum(Y(:,desno))/nr<200)
    desno=desno-1;
end
X=X(poc:kraj,levo:desno);
% figure(2);imshow(X)
Y=Y(poc:kraj,levo:desno);
[nr,nc]=size(X);

% isecanje belih segmenata
poc=1;
while (poc<nr) && (sum(Y(poc,:))/nc>250) || ...
        ((sum(Y(poc,:))/nc<250) && (sum(Y(poc+1,:))/nc>250))
    poc=poc+1;
end

kraj=nr;
while (kraj>1) && (sum(Y(kraj,:))/nc>250) || ...
        ((sum(Y(kraj,:))/nc<250) && (sum(Y(kraj-1,:))/nc>250))
    kraj=kraj-1;
end

levo=1;
while (levo<nc) && (sum(Y(:,levo))/nr>250) || ...
        ((sum(Y(:,levo))/nr<250) && (sum(Y(:,levo+1))/nr>250))
    levo=levo+1;
end

desno=nc;
while (desno>1) && (sum(Y(:,desno))/nr>250) || ...
        ((sum(Y(:,desno))/nr<250) && (sum(Y(:,desno-1))/nr>250))
    desno=desno-1;
end

X=X(poc:kraj,levo:desno);
Y=Y(poc:kraj,levo:desno);
if a==1
    figure(50);imshow(Y/255); pause
end
[nr,nc]=size(X);

% presek po 1/3 
lvl1 = round(nr/3); % 1/3
razmak1 = zeros(10,1); % ako je razmak izmedju prethodne i sledece linije mali, 
noUp1 = 1; % koliko je puta presecena crna linija na novou1
start1 = 0; % fleg - da znamo da li smo poceli sa crnom linijom
for i = 1:nc
    if ((Y(lvl1,i)==255) && (start1 == 0))
        start1 = 1;
        noUp1 = noUp1 + 1;
    end
    if ((Y(lvl1,i)==0) && (start1 == 0))
        razmak1(noUp1) = razmak1(noUp1) + 1;
    end
    
    if ((Y(lvl1,i)==0) && (start1 == 1))
        start1 = 0; 
        razmak1(noUp1) = 0;
    end 
end

razmak1(razmak1<4) = 0;
razmak1(razmak1>=4) = 1;
noUp1 = sum(razmak1)+1;

% mean u sredini slike
P(1,1) = mean(mean(Y(round(1/3*nr):round(2/3*nr),...
    round(1/3*nc):round(2/3*nc))))/255;
% mean u sredisnjoj donjoj cetvrtini
P(2,1) = mean(mean(Y((round(3/4*nr):nr),round(1/3*nc):round(2/3*nc))))/255;
% odnos gornje i donje polovine
P(3,1) = sum(sum(Y(1:round(1/2*nr),:)))/sum(sum(Y((round(1/2*nr)+1):nr,:)));
% odnos leve i desne 1/3
P(4,1) = mean(mean(Y(:,1:round(1/3*nc))))/mean(mean(Y(:,(round(2/3*nc)+1):nc)));
% broj preseka na 1/3 slike (odozgo gledano)
P(5,1) = noUp1;


end
