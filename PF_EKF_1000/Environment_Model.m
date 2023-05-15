%ENVIRONMENT MODEL
%Sviluppato dal gruppo C4 del Team C
% Membri :
% - Mucci Beltrami Marco
% - Rocco Luigi Gaspare Burdo
% - Giacomo Bertocchini
% - Giulio Pisaneschi

% Questo programma plotta una traiettoria ideale del veicolo, partendo
% dal punto iniziale in superficie si porta al punto di inizio
% pattugliamento in superficie, si immerge, segue poi il percorso indicato
% dal file missionC. 
% Oltre questo plotta anche il percorso realmente fatto dal veicolo. 
% L'implementazione dei sonar è fatta tramite una loro schematizzazione
% a 'cono' prendendo 9 punti come i punti centrali ed i vertici di un rettangolo.

close all
%% Inizializzazione
% In questa sezione di codice viene inizializzato il programma, riportando i
% punti dei corner ed il punto di inizio missione in una terna NED centrata 
% nello spigolo C. Viene inoltre inizializzata la finestra dove avverrà il
% plot.

wgs84 = wgs84Ellipsoid;
[xA,yA,zA] = geodetic2ned(cornerA(1),cornerA(2),0,cornerC(1),cornerC(2),0,wgs84);
[xB,yB,zB] = geodetic2ned(cornerB(1),cornerB(2),0,cornerC(1),cornerC(2),0,wgs84);
[xC,yC,zC] = geodetic2ned(cornerC(1),cornerC(2),0,cornerC(1),cornerC(2),0,wgs84);
[xD,yD,zD] = geodetic2ned(cornerD(1),cornerD(2),0,cornerC(1),cornerC(2),0,wgs84);
[xiP,yiP,ziP]=geodetic2ned(initPoint(1),initPoint(2),-0.485,cornerC(1),cornerC(2),0,wgs84);
ziP=-ziP;
iPNed=[xiP,yiP,ziP];
Env= figure('Name','Team C: Environment Model','NumberTitle','off');

%% Plot del bacino
%In questa sezione viene plottato il bacino.L'altezza è variabile in
%funzione della profondità fisata nel file di inizializzazione.

%creazione dei vettori utili al plot delle pareti del bacino
x1 = [xA;xB]; 
y1 = [yA;yB];
x2 = [xB;xC];
y2 = [yB;yC];
x3 = [xC;xD];
y3 = [yC;yD];
x4 = [xD;xA];
y4 = [yD;yA];
x = [x1 x2 x3 x4];
y = [y1 y2 y3 y4];

if depth>25 %Profondità bacino
    h1=-depth-5;
else
h1 = -30; 
end

z=[0;0;h1;h1];
Color1=[222,184,135]; %Impostazione del colore del bacino in RGB
Color=Color1/256;
h = fill3(  [yA yB yB yA],[xA xB xB xA],z,Color,...
            [yB yC yC yB],[xB xC xC xB],z,Color,...
            [yC yD yD yC],[xC xD xD xC],z,Color,...
            [yD yA yA yD],[xD xA xA xD],z,Color);
        
set(h,'FaceAlpha',0.7) ;

%plot del piano dell'acqua
s=patch([yA yB yC yD],[xA xB xC xD], [0 0 0 0],'FaceAlpha',0.5);
s.EdgeColor = 'none';

% %Uncomment se in NED
axis equal
xlabel('East');
ylabel('North');
zlabel('Down');

text(yC,xC+8,'C');
text(yD,xD+8,'D');
text(yA,xA+8,'A');
text(yB,xB+8,'B');



hold on; %plot del punto di inizio ispezione
plot3(yiP, xiP,ziP,'ro', 'MarkerSize', 10);
hold on

%% Calcolo dei punti notevoli
% I punti notevoli sono i punti sugli spigoli e 
% i punti di inizio ispezione sulla parete. 
% Questi si trovano a una distanza dal muro data dal file di inizializzazione
% e fungeranno da vertici dei segmenti ideali dell’ispezione.

%Legenda:pip_ è il punto di inizio pattugliamento, pip_prof è il punto di inizio
%pattugliamento in profondità

%calcolo lunghezze dei lati
AB=sqrt((xA-xB)^2+(yA-yB)^2); %lunghezza AB
BC=sqrt((xB-xC)^2+(yB-yC)^2); %lunghezza BC
CD=sqrt((xC-xD)^2+(yC-yD)^2); %lunghezza CD
DA=sqrt((xD-xA)^2+(yD-yA)^2); %lughezza DA

dW=distanceFromWall;
%Legenda:pc_ è il punto centrale di ogni parete
pcAB=[(xB+xA)/2,(yB+yA)/2]; 
pcBC=[(xB+xC)/2,(yB+yC)/2];
pcCD=[(xC+xD)/2,(yC+yD)/2];
pcDA=[(xD+xA)/2,(yD+yA)/2];

% Calcolo angoli [rad] fra le pareti

%In questi calcoli viene sfruttata la formula per i triangoli rettangoli
%BC*cos(thetaBC)=yB, si ricorda che le y sono le ascisse. 
thetaBC= acos(yB/BC);
%thetaBCdeg=rad2deg(thetaBC);

%CD*cos(thetaCD)=yD
thetaCD=acos(yD/CD);
%thetaCDdeg=rad2deg(thetaCD);

%AB*cos(gammaAB)=(yB-yA)
gammaAB=acos((yB-yA)/AB);
thetaAB=pi-gammaAB;
%thetaABdeg=rad2deg(thetaAB);

%DA*cos(gammaDA)=(yA-yD)
gammaDA=acos((yA-yD)/DA);
%gammaDAdeg=rad2deg(gammaDA);

% questa sequenza di codice genera il punto di inizio pattugliamento a
% distanza dW dalla parete AB
dbpc=sqrt((pcAB(1)-xB)^2+(pcAB(2)-yB)^2); %distanza dal punto centrale di AB a B
hpb=sqrt(dbpc^2 + dW^2); %ipotenusa B-puntocentrale-DW
alfa1=-acos(dbpc/hpb);%angolo fra hp e lato B-PC
gammaab=pi-(gammaAB+alfa1);
pcab1=[hpb*sin(gammaab),hpb*cos(gammaab)];
pipAB=[xB+pcab1(1),yB+pcab1(2),0];
pipABprof=[pipAB(1),pipAB(2),-depth];

%questa sequenza di codice genera il punto di inizio pattugliamento a
% distanza dW dalla parete BC
dcpc=sqrt(pcBC(1)^2+pcBC(2)^2);
hpc=sqrt(dcpc^2 + dW^2);
alfabc=acos(dcpc/hpc);
gammaBC=thetaBC+alfabc;
pipBC=[hpc*sin(gammaBC),hpc*cos(gammaBC),0];
pipBCprof=[pipBC(1),pipBC(2),-depth];

%questa sequenza di codice genera il punto di inizio pattugliamento a
% distanza dW dalla parete CD
dcdpc=sqrt(pcCD(1)^2+pcCD(2)^2);
hpcd=sqrt(dcdpc^2 + dW^2);
alfacd=acos(dcdpc/hpcd);
gammaCD=thetaCD-alfacd;
pipCD=[hpcd*sin(gammaCD),hpcd*cos(gammaCD),0];
pipCDprof=[pipCD(1),pipCD(2),-depth];

%questa sequenza di codice genera il punto di inizio pattugliamento a
% distanza dW dalla parete DA
dapc=sqrt((pcDA(1)-xA)^2+(pcDA(2)-yA)^2);
hpda=sqrt(dapc^2 + dW^2);
alfada=acos(dapc/hpda);
gammada=pi+(gammaDA+alfada);
pcda1=[hpda*sin(gammada),hpda*cos(gammada)];
pipDA=[xA+pcda1(1),yA+pcda1(2),0];
pipDAprof=[pipDA(1),pipDA(2),-depth];

% Questa sequenza di codice genera i punti sugli spigoli, quei punti cioè
% che si trovano a distanza dW da ciascuna parete

% calcolo punto sullo spigolo B
BC1=BC-dW;
hypspigoloB=sqrt(BC1^2+dW^2);
beta1=acos(BC1/hypspigoloB);
betaB=thetaBC+beta1;
spigB=[hypspigoloB*sin(betaB),hypspigoloB*cos(betaB),-depth];

%calcolo punto sullo spigolo C
hypspigoloC=sqrt(dW^2+dW^2);
beta2=acos(dW/hypspigoloC);
betaC=thetaBC+beta2;
spigC=[hypspigoloC*sin(betaC),hypspigoloC*cos(betaC),-depth];

%calcolo punto sullo spigolo D
CD1=CD-dW;
hypspigoloD=sqrt(CD1^2+dW^2);
beta3=-acos(CD1/hypspigoloD);
betaD=thetaCD+beta3;
spigD=[hypspigoloD*sin(betaD),hypspigoloD*cos(betaD),-depth];

%calcolo punto sullo spigolo A
AB1=AB-dW;
hypspigoloA=sqrt(AB1^2+dW^2);
beta4=acos(AB1/hypspigoloA);
betaA=thetaAB+beta4;
spigA1=[hypspigoloA*sin(betaA),hypspigoloA*cos(betaA)];
spigA=[spigA1(1)+xB,spigA1(2)+yB,-depth];

%In questa sequenza vengono plottati i punti appena calcolati

plot3(pipCD(2),pipCD(1),-depth,'c*', 'MarkerSize', 10);
hold on
plot3(pipAB(2),pipAB(1),-depth,'c*', 'MarkerSize', 10);
hold on
plot3(pipBC(2),pipBC(1),-depth,'c*', 'MarkerSize', 10);
hold on
plot3(pipBCprof(2),pipBCprof(1),-depth,'c*', 'MarkerSize', 10);
hold on
plot3(pipDA(2),pipDA(1),-depth,'c*', 'MarkerSize', 10);
hold on
plot3(spigB(2),spigB(1),spigB(3),'c*', 'MarkerSize', 10);
hold on
plot3(spigC(2),spigC(1),spigC(3),'c*', 'MarkerSize', 10);
hold on
plot3(spigD(2),spigD(1),spigD(3),'c*', 'MarkerSize', 10);
hold on
plot3(spigA(2),spigA(1),spigA(3),'c*', 'MarkerSize', 10);

%% Calcolo traiettoria ideale
%In questa sequenza di codice vengono generati i vettori lungo i quali si
%svolgera' il plot ideale. Una volta generati verranno riordinati in
%funzione dei parametri di missione.

%Legenda: ip=Initial Point
hold on
Nc=5; %numero di campioni per tratto di ispezione

%def parte iniziale ip->AB
xAB=linspace(iPNed(1),pipAB(1),Nc);
yAB=linspace(iPNed(2),pipAB(2),Nc);
zAB=linspace(iPNed(3),pipAB(3),Nc);
iPpipAB=[xAB;yAB;zAB];

xABp=linspace(pipAB(1),pipABprof(1),Nc);
yABp=linspace(pipAB(2),pipABprof(2),Nc);
zABp=linspace(pipAB(3),pipABprof(3),Nc);
pipabprof=[xABp;yABp;zABp];

%def parte iniziale ip->BC
xBC=linspace(iPNed(1),pipBC(1),Nc);
yBC=linspace(iPNed(2),pipBC(2),Nc);
zBC=linspace(iPNed(3),pipBC(3),Nc);
iPpipBC=[xBC;yBC;zBC];

xBCp=linspace(pipBC(1),pipBCprof(1),Nc);
yBCp=linspace(pipBC(2),pipBCprof(2),Nc);
zBCp=linspace(pipBC(3),pipBCprof(3),Nc);
pipbcprof=[xBCp;yBCp;zBCp];

%def parte iniziale ip->CD
xCD=linspace(iPNed(1),pipCD(1),Nc);
yCD=linspace(iPNed(2),pipCD(2),Nc);
zCD=linspace(iPNed(3),pipCD(3),Nc);
iPpipCD=[xCD;yCD;zCD];

xCDp=linspace(pipCD(1),pipCDprof(1),Nc);
yCDp=linspace(pipCD(2),pipCDprof(2),Nc);
zCDp=linspace(pipCD(3),pipCDprof(3),Nc);
pipcdprof=[xCDp;yCDp;zCDp];

%def parte iniziale ip->DA
xDA=linspace(iPNed(1),pipDA(1),Nc);
yDA=linspace(iPNed(2),pipDA(2),Nc);
zDA=linspace(iPNed(3),pipDA(3),Nc);
iPpipDA=[xDA;yDA;zDA];

xDAp=linspace(pipDA(1),pipDAprof(1),Nc);
yDAp=linspace(pipDA(2),pipDAprof(2),Nc);
zDAp=linspace(pipDA(3),pipDAprof(3),Nc);
pipdaprof=[xDAp;yDAp;zDAp];

% riordinamento dei vettori per il plot in funzione dei parametri
% della missione

Vect=[spigA;pipABprof;spigB;pipBCprof;spigC;pipCDprof;spigD;pipDAprof];

l=length(Vect);
Vectsorted=zeros(l,3);

if(initInspectedSide==1)
    start=2;
elseif(initInspectedSide==2)
    start=4;
elseif(initInspectedSide==3)
    start=6;
elseif(initInspectedSide==4)
    start=8;
end

if  strcmpi(inspectionDirection,'clockwise')
for i=0:(l-1)
    if (start+i)>8
    Vectsorted(i+1,:)=Vect(abs(l-(i+start)),:);
    else
        Vectsorted(i+1,:)=Vect(start+i,:);
    end
end
elseif strcmpi(inspectionDirection,'counterclockwise')
    for i=0:(l-1)
        if (start-i)<=0
            Vectsorted(i+1,:)=Vect(l+start-i,:);
        else
            Vectsorted(i+1,:)=Vect(start-i,:);
        end
    end
end

if(initInspectedSide==1)
    xtraj1=iPpipAB(1,:);
    xtraj2=pipabprof(1,:);
    ytraj1=iPpipAB(2,:);
    ytraj2=pipabprof(2,:);
    ztraj1=iPpipAB(3,:);
    ztraj2=pipabprof(3,:);
elseif(initInspectedSide==2)
    xtraj1=iPpipBC(1,:);
    xtraj2=pipbcprof(1,:);
    ytraj1=iPpipBC(2,:);
    ytraj2=pipbcprof(2,:);
    ztraj1=iPpipBC(3,:);
    ztraj2=pipbcprof(3,:);
elseif(initInspectedSide==3)
    xtraj1=iPpipCD(1,:);
    xtraj2=pipcdprof(1,:);
    ytraj1=iPpipCD(2,:);
    ytraj2=pipcdprof(2,:);
    ztraj1=iPpipCD(3,:);
    ztraj2=pipcdprof(3,:);
elseif(initInspectedSide==4)
    xtraj1=iPpipDA(1,:);
    xtraj2=pipdaprof(1,:);
    ytraj1=iPpipDA(2,:);
    ytraj2=pipdaprof(2,:);
    ztraj1=iPpipDA(3,:);
    ztraj2=pipdaprof(3,:);
end


% costruzione dei vettori per il plot della traiettoria ideale

xtraj=[xtraj1',xtraj2',linspace(Vectsorted(1),Vectsorted(2),Nc)',linspace(Vectsorted(2),Vectsorted(3),Nc)',linspace(Vectsorted(3),Vectsorted(4),Nc)',...
    linspace(Vectsorted(4),Vectsorted(5),Nc)',linspace(Vectsorted(5),Vectsorted(6),Nc)',linspace(Vectsorted(6),Vectsorted(7),Nc)',...
    linspace(Vectsorted(7),Vectsorted(8),Nc)',linspace(Vectsorted(8),Vectsorted(1),Nc)'];
ytraj=[ytraj1',ytraj2',linspace(Vectsorted(9),Vectsorted(10),Nc)',linspace(Vectsorted(10),Vectsorted(11),Nc)',linspace(Vectsorted(11),Vectsorted(12),Nc)',...
    linspace(Vectsorted(12),Vectsorted(13),Nc)',linspace(Vectsorted(13),Vectsorted(14),Nc)',linspace(Vectsorted(14),Vectsorted(15),Nc)',...
    linspace(Vectsorted(15),Vectsorted(16),Nc)',linspace(Vectsorted(16),Vectsorted(9),Nc)'];
ztraj=[ztraj1',ztraj2',linspace(Vectsorted(17),Vectsorted(18),Nc)',linspace(Vectsorted(18),Vectsorted(19),Nc)',linspace(Vectsorted(19),Vectsorted(20),Nc)',...
    linspace(Vectsorted(20),Vectsorted(21),Nc)',linspace(Vectsorted(21),Vectsorted(22),Nc)',linspace(Vectsorted(22),Vectsorted(23),Nc)',...
    linspace(Vectsorted(23),Vectsorted(24),Nc)',linspace(Vectsorted(24),Vectsorted(17),Nc)'];

len1=size(xtraj);

%Plot della traiettoria ideale

for j=1:(len1(2))
    for i=1:(len1(1)-1)
    h=plot3(ytraj(i:i+1,j),xtraj(i:i+1,j),ztraj(i:i+1,j),'r','linewidth',1);
    pause(0.05) 
    end
end
%% Calcolo dei piani delle pareti e dell'acqua
% Questa sezione ricava i coefficienti dei piani delle pareti.
% L'equazione del piano corrisponde a quella della retta perché il piano
% delle pareti è verticale, per cui si risolve un sistema lineare:
% l'equazione del piano deve soddisfare il passaggio da due punti.
% Il parametro a è libero, lo poniamo a=1, il parametro c relativo all'asse
% z sarà sempre pari a zero in quanto le pareti non sono inclinate.

% Questo è il sistema da risolvere per trovare i parametri delle rette
% passanti per i punti
m_AB = - (yB- yA) /( xB - xA );   %AB
q_AB = - yA - xA*m_AB;

m_BC = - (yC-yB)/(xC-xB);         %BC
q_BC = - yB - xB*m_BC;

m_CD = - (yD - yC) / (xD - xC);   %CD
q_CD = - yC - xC*m_CD;

m_DA = - (yA - yD)/(xA - xD);     %DA
q_DA = -yD - xD*m_DA;

% Noti i parametri precedenti, possiamo ricavare i parametri dei piani
a_AB = 1;
b_AB = double(m_AB);
c_AB = 0;
d_AB = double(q_AB);
 

a_BC = 1;
b_BC = double(m_BC);
c_BC = 0;
d_BC = double(q_BC);

a_CD = 1;
b_CD = double(m_CD);
c_CD = 0;
d_CD = double(q_CD);

a_DA = 1;
b_DA = double(m_DA);
c_DA = 0;
d_DA = double(q_DA);

% Parametri dei piani delle pareti
plane_AB = [a_AB b_AB c_AB d_AB];
plane_BC = [a_BC b_BC c_BC d_BC];
plane_CD = [a_CD b_CD c_CD d_CD];
plane_DA = [a_DA b_DA c_DA d_DA];
% Parametri del piano dell'acqua
plane_W =[0 0 1 0];
hold on
%% Implementazione dei Sonar

% Per ogni punto del pattugliamento e
% per la sua orientazione vengono calcolate le intersezioni del fascio del
% sonar con le pareti e le relative distanze. Questi valori vengono
% salvati in matrici apposite che verranno riutilizzate per il plot.

% Questa sequenza di codice salva gli output del simulink, verrà utilizzato
% per il plot della traiettoria realmente seguita.
xreal=out.Lon_ts.Data;
yreal=out.Lat_ts.Data;
zreal=-out.Depth_ts.Data;
Ncr=length(xreal);

% Le variabili p_ conterranno al proprio interno le coordinate delle
% intersezioni fra le pareti e la retta del sonar. La loro dimensione
% deriva dal fatto che per ogni istante si avranno 3 coordinate (x,y,z), e
% per ogni intercetta verrà simulata la presenza di un 'cono' composto da 9
% punti: il punto centrale, i punti sugli spigoli di un rettangolo con lati
% uguali al fascio del sonar ed i punti al centro di tali lati. 
% Una sua rappresentazione 2d è la seguente:
% x   x   x
% x   x   x
% x   x   x

cono=zeros(3,3); %cono conterrà i 9 punti del sonar
p1=zeros(Ncr,3,3,3);
p2=zeros(Ncr,3,3,3);
p3=zeros(Ncr,3,3,3);
Ns=3;

WBr=deg2rad(sens_BeamWidth/2); %larghezza del semifascio del cono

distsonarprua=zeros(1,Ncr);
distsonardx=zeros(1,Ncr);
distsonarsx=zeros(1,Ncr);

% In questa simulazione c'è una precisazione da fare, per costruzione del
% plot di matlab, gli assi x ed y sono stati invertiti. Il Bacino è stato
% cioè plottato con yEast sulle ascisse e xNorth sulle ordinate. In
% conseguenza a questo fatto i versori dei sonar saranno :
% - Versore del Sonar di Prua = [0 1 0]
% - Versore del Sonar di Destra = [1 0 0]
% - Versore del Sonar di Sinistra = [-1 0 0]


% Il codice seguente ripete lo stesso procedimento per i 3 sonar,
% sostanzialmente si ha il calcolo delle intercette di ciascun ramo del 
% 'cono', quando queste saranno sul piano dell'acqua verrà restuito 'NaN',
% se invece sono su una delle pareti del bacino verrà
% restituito il valore della sua distanza. Nel caso in cui si abbiano
% almeno 5 NaN (su 9 valori) verrà restituito NaN, altrimenti verrà
% calcolata la media escludendo i valori non validi. Questo valore della
% distanza viene salvato nella variabile 'distsonar_'

versore_son_body = [0 1 0]'; % Implementazione cono nel Sonar di Prua
for k=1:Ncr
    Roll=out.Roll_ts.Data(k);
    Pitch=out.Pitch_ts.Data(k);
    Yaw=out.Yaw_ts.Data(k);
    
% creazione del cono 
R=Roll;
P=linspace(Pitch-WBr,Pitch+WBr,Ns);
Y=linspace(Yaw-WBr, Yaw+WBr,Ns);

 for j=1:Ns
  for i=1:Ns
               
orient= [P(j) R Y(i)]'; 

% Calcolo del Jacobiano
R_x = [1        0              0;
      0 cos(orient(1))  -sin(orient(1));
      0 sin(orient(1)) cos(orient(1))];
  
R_y = [cos(orient(2)) 0 sin(orient(2));
          0           1         0;
       -sin(orient(2)) 0 cos(orient(2))];
   
R_z = [cos(orient(3)) -sin(orient(3)) 0;
      sin(orient(3)) cos(orient(3)) 0;
          0                 0        1];
jaco = (R_z'*R_y');
jacob = jaco*R_x';
sonar_son_NED = jacob*versore_son_body;

    pos=[xreal(k,1),yreal(k,1),zreal(k,1)];
    
    %Calcolo delle intercette con i piani
    
    t_AB = -(pos(1)+plane_AB(2)*pos(2)+plane_AB(4))/(plane_AB(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_BC = -(pos(1)+plane_BC(2)*pos(2)+plane_BC(4))/(plane_BC(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_CD = -(pos(1)+plane_CD(2)*pos(2)+plane_CD(4))/(plane_CD(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_DA = -(pos(1)+plane_DA(2)*pos(2)+plane_DA(4))/(plane_DA(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_W = -(plane_W(3)*pos(3)/(plane_W(3)*sonar_son_NED(3)+sonar_son_NED(3)));
    % Coordinate delle intersezioni con i piani
    intx_AB(1) = pos(1) + sonar_son_NED(1)*double(t_AB);
    intx_AB(2) = pos(2) + sonar_son_NED(2)*double(t_AB);
    intx_AB(3) = pos(3) + sonar_son_NED(3)*double(t_AB);
    
    intx_BC(1) = pos(1) + sonar_son_NED(1)*double(t_BC);
    intx_BC(2) = pos(2) + sonar_son_NED(2)*double(t_BC);
    intx_BC(3) = pos(3) + sonar_son_NED(3)*double(t_BC);
    
    intx_CD(1) = pos(1) + sonar_son_NED(1)*double(t_CD);
    intx_CD(2) = pos(2) + sonar_son_NED(2)*double(t_CD);
    intx_CD(3) = pos(3) + sonar_son_NED(3)*double(t_CD);
    
    intx_DA(1) = pos(1) + sonar_son_NED(1)*double(t_DA);
    intx_DA(2) = pos(2) + sonar_son_NED(2)*double(t_DA);
    intx_DA(3) = pos(3) + sonar_son_NED(3)*double(t_DA);
    
    intx_W(1) = pos(1) + sonar_son_NED(1)*double(t_W);
    intx_W(2) = pos(2) + sonar_son_NED(2)*double(t_W);
    intx_W(3) = pos(3) + sonar_son_NED(3)*double(t_W);
    
    % I valori negativi vengono considerati non validi
    
    if t_AB < 0
        t_AB = inf;
    end
    if t_BC < 0
        t_BC = inf;
    end
    if t_CD < 0
        t_CD = inf;
    end
    if t_DA < 0
        t_DA = inf;
    end
    if t_W < 0
        t_W = inf;
    end
    
    % Si va a considerare corretto il valore minore, sarà il primo piano ad
    % essere intercettato dalla retta del sonar
    dist = min([t_AB t_BC t_CD t_DA t_W]);
    
hold on

%salvo in p1 le coordinate dell'intercetta
if dist == t_AB
p1(k,:,i,j)= intx_AB;

elseif dist == t_BC
p1(k,:,i,j)=intx_BC;

elseif dist == t_CD
p1(k,:,i,j)=intx_CD;

elseif dist == t_DA
p1(k,:,i,j)=intx_DA;

elseif dist == t_W
p1(k,:,i,j)=intx_W;
end

    if dist==t_W % se incontro l'acqua metto NaN nella misura delle distanze
     cono(i,j)=NaN;
    else
     cono(i,j)=dist;
    end
           
   end
  end
    
    contr=isnan(cono); %controlla la presenza di NaN nel vettore di misure
    sum1=sum(contr,'all'); % conta il numero di NaN presenti
    
    if sum1>4 %incontro un numero di Nan superiore a 4 considero la misura non valida
        distsonarprua(k)=NaN;
    else
        distsonarprua(k)=nanmean(cono,'all');
    end    
    
end

versore_son_body = [1 0 0]'; % Implementazione cono nel Sonar di Destra

for k=1:Ncr
 
            
    Roll=out.Roll_ts.Data(k);
    Pitch=out.Pitch_ts.Data(k);
    Yaw=out.Yaw_ts.Data(k);
    
R=linspace(Roll-WBr,Roll+WBr,Ns);
P=Pitch;
Y=linspace(Yaw-WBr, Yaw+WBr,Ns);

     for j=1:Ns
        for i=1:Ns
orient= [P R(j) Y(i)]'; 

R_x = [1        0              0;
      0 cos(orient(1))  -sin(orient(1));
      0 sin(orient(1)) cos(orient(1))];
  
R_y = [cos(orient(2)) 0 sin(orient(2));
          0           1         0;
       -sin(orient(2)) 0 cos(orient(2))];
   
R_z = [cos(orient(3)) -sin(orient(3)) 0;
      sin(orient(3)) cos(orient(3)) 0;
          0                 0        1];
jaco = (R_z'*R_y');
jacob = jaco*R_x';
sonar_son_NED = jacob*versore_son_body;

    pos=[xreal(k,1),yreal(k,1),zreal(k,1)];
    
    t_AB = -(pos(1)+plane_AB(2)*pos(2)+plane_AB(4))/(plane_AB(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_BC = -(pos(1)+plane_BC(2)*pos(2)+plane_BC(4))/(plane_BC(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_CD = -(pos(1)+plane_CD(2)*pos(2)+plane_CD(4))/(plane_CD(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_DA = -(pos(1)+plane_DA(2)*pos(2)+plane_DA(4))/(plane_DA(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_W = -(plane_W(3)*pos(3)/(plane_W(3)*sonar_son_NED(3)+sonar_son_NED(3)));
    
    intx_AB(1) = pos(1) + sonar_son_NED(1)*double(t_AB);
    intx_AB(2) = pos(2) + sonar_son_NED(2)*double(t_AB);
    intx_AB(3) = pos(3) + sonar_son_NED(3)*double(t_AB);
    
    intx_BC(1) = pos(1) + sonar_son_NED(1)*double(t_BC);
    intx_BC(2) = pos(2) + sonar_son_NED(2)*double(t_BC);
    intx_BC(3) = pos(3) + sonar_son_NED(3)*double(t_BC);
    
    intx_CD(1) = pos(1) + sonar_son_NED(1)*double(t_CD);
    intx_CD(2) = pos(2) + sonar_son_NED(2)*double(t_CD);
    intx_CD(3) = pos(3) + sonar_son_NED(3)*double(t_CD);
    
    intx_DA(1) = pos(1) + sonar_son_NED(1)*double(t_DA);
    intx_DA(2) = pos(2) + sonar_son_NED(2)*double(t_DA);
    intx_DA(3) = pos(3) + sonar_son_NED(3)*double(t_DA);
    
    intx_W(1) = pos(1) + sonar_son_NED(1)*double(t_W);
    intx_W(2) = pos(2) + sonar_son_NED(2)*double(t_W);
    intx_W(3) = pos(3) + sonar_son_NED(3)*double(t_W);
    
    if t_AB < 0
        t_AB = inf;
    end
    if t_BC < 0
        t_BC = inf;
    end
    if t_CD < 0
        t_CD = inf;
    end
    if t_DA < 0
        t_DA = inf;
    end
    if t_W < 0
        t_W = inf;
    end
           
    dist = min([t_AB t_BC t_CD t_DA t_W]);

hold on

if dist == t_AB
p2(k,:,i,j)= intx_AB;

elseif dist == t_BC
p2(k,:,i,j)=intx_BC;

elseif dist == t_CD
p2(k,:,i,j)=intx_CD;

elseif dist == t_DA
p2(k,:,i,j)=intx_DA;

elseif dist == t_W
p2(k,:,i,j)=intx_W;
end

           if dist==t_W %se incontro l'acqua metto NaN nella misura delle distanze
               cono(i,j)=NaN;
           else
               cono(i,j)=dist;
           end
        end
     end
    
    contr=isnan(cono); %controlla la presenza di NaN nel vettore di misure
    sum1=sum(sum(contr)); % conta il numero di NaN presenti
    
    if sum1>4 %incontro un numero di Nan superiore a 4
        distsonardx(1,k)=NaN;
    else
        distsonardx(1,k)=nanmean(cono,'all');
    end    
    
end

versore_son_body = [-1 0 0]'; % Implementazione cono nel Sonar di Sinistra
for k=1:Ncr
     
            
    Roll=out.Roll_ts.Data(k);
    Pitch=out.Pitch_ts.Data(k);
    Yaw=out.Yaw_ts.Data(k);
    
R=linspace(Roll-WBr,Roll+WBr,Ns);
P=Pitch;
Y=linspace(Yaw-WBr, Yaw+WBr,Ns);

     for i=1:Ns
        for j=1:Ns
orient= [P R(j) Y(i)]'; 

R_x = [1        0              0;
      0 cos(orient(1))  -sin(orient(1));
      0 sin(orient(1)) cos(orient(1))];
  
R_y = [cos(orient(2)) 0 sin(orient(2));
          0           1         0;
       -sin(orient(2)) 0 cos(orient(2))];
   
R_z = [cos(orient(3)) -sin(orient(3)) 0;
      sin(orient(3)) cos(orient(3)) 0;
          0                 0        1];
jaco = (R_z'*R_y');
jacob = jaco*R_x';
sonar_son_NED = jacob*versore_son_body;

    pos=[xreal(k,1),yreal(k,1),zreal(k,1)];
    
    t_AB = -(pos(1)+plane_AB(2)*pos(2)+plane_AB(4))/(plane_AB(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_BC = -(pos(1)+plane_BC(2)*pos(2)+plane_BC(4))/(plane_BC(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_CD = -(pos(1)+plane_CD(2)*pos(2)+plane_CD(4))/(plane_CD(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_DA = -(pos(1)+plane_DA(2)*pos(2)+plane_DA(4))/(plane_DA(2)*sonar_son_NED(2)+sonar_son_NED(1));
    t_W = -(plane_W(3)*pos(3)/(plane_W(3)*sonar_son_NED(3)+sonar_son_NED(3)));
    
    intx_AB(1) = pos(1) + sonar_son_NED(1)*double(t_AB);
    intx_AB(2) = pos(2) + sonar_son_NED(2)*double(t_AB);
    intx_AB(3) = pos(3) + sonar_son_NED(3)*double(t_AB);
    
    intx_BC(1) = pos(1) + sonar_son_NED(1)*double(t_BC);
    intx_BC(2) = pos(2) + sonar_son_NED(2)*double(t_BC);
    intx_BC(3) = pos(3) + sonar_son_NED(3)*double(t_BC);
    
    intx_CD(1) = pos(1) + sonar_son_NED(1)*double(t_CD);
    intx_CD(2) = pos(2) + sonar_son_NED(2)*double(t_CD);
    intx_CD(3) = pos(3) + sonar_son_NED(3)*double(t_CD);
    
    intx_DA(1) = pos(1) + sonar_son_NED(1)*double(t_DA);
    intx_DA(2) = pos(2) + sonar_son_NED(2)*double(t_DA);
    intx_DA(3) = pos(3) + sonar_son_NED(3)*double(t_DA);
    
    intx_W(1) = pos(1) + sonar_son_NED(1)*double(t_W);
    intx_W(2) = pos(2) + sonar_son_NED(2)*double(t_W);
    intx_W(3) = pos(3) + sonar_son_NED(3)*double(t_W);
    
    if t_AB < 0
        t_AB = inf;
    end
    if t_BC < 0
        t_BC = inf;
    end
    if t_CD < 0
        t_CD = inf;
    end
    if t_DA < 0
        t_DA = inf;
    end
    if t_W < 0
        t_W = inf;
    end
           
    dist = min([t_AB t_BC t_CD t_DA t_W]);
    
hold on

if dist == t_AB
p3(k,:,i,j)= intx_AB;

elseif dist == t_BC
p3(k,:,i,j)=intx_BC;

elseif dist == t_CD
p3(k,:,i,j)=intx_CD;

elseif dist == t_DA
p3(k,:,i,j)=intx_DA;

elseif dist == t_W
p3(k,:,i,j)=intx_W;
end
           if dist==t_W %se incontro l'acqua metto NaN nella misura delle distanze
               cono(i,j)=NaN;
           else
               cono(i,j)=dist;
           end

        end
     end
        contr=isnan(cono); %controlla la presenza di NaN nel vettore di misure
    sum1=sum(sum(contr)); % conta il numero di NaN presenti
    
    if sum1>4 %incontro un numero di Nan superiore a 4
        distsonarsx(k)=NaN;
    else
        distsonarsx(k)=nanmean(cono,'all');
    end   
end

%% Plot Complessivo
% In questo ultimo ciclo for vengono plottati :
% - La traiettoria seguita dal veicolo (linea Blu)
% - Il punto in cui si trova il veicolo (punto Viola)
% - i coni dei sonar
% - i punti di intersezione con i piani
for i=1:(Ncr-1)
     pos=[xreal(i,1),yreal(i,1),zreal(i,1)];
     % plot traiettoria e punto
      plot3(xreal(i:i+1),yreal(i:i+1),zreal(i:i+1),'b','linewidth',2)
      deep_purple=plot3(xreal(i),yreal(i),zreal(i),'o','MarkerFaceColor',[0.5,0,0.5],'MarkerSize',10);
      
      %plot punto centrale del sonar di prua (giallo)
      c5t=plot3([pos(1) p1(i,1,2,2)],[pos(2) p1(i,2,2,2)],[pos(3) p1(i,3,2,2)],'y'); %prua
      c5p1=plot3(p1(i,1,2,2),p1(i,2,2,2),p1(i,3,2,2),'y*');
      t1=text(p1(i,1,2,2)+1,p1(i,2,2,2)+1,num2str(distsonarprua(i)));
      
      %plot punto centrale del sonar di destra (verde)
      d5t=plot3([pos(1) p2(i,1,2,2)],[pos(2) p2(i,2,2,2)],[pos(3) p2(i,3,2,2)],'g'); %dx
      d5p=plot3(p2(i,1,2,2),p2(i,2,2,2),p2(i,3,2,2),'g*');
      t2=text(p2(i,1,2,2)+1,p2(i,2,2,2)+1,num2str(distsonardx(i)));
      
      %plot punto centrale del sonar di sinistra (magenta)
      s5t=plot3([pos(1) p3(i,1,2,2)],[pos(2) p3(i,2,2,2)],[pos(3) p3(i,3,2,2)],'m'); %sx
      s5p=plot3(p3(i,1,2,2),p3(i,2,2,2),p3(i,3,2,2),'m*');
      t3=text(p3(i,1,2,2)+1,p3(i,2,2,2)+1,num2str(distsonarsx(i)));
      
    %Commentare la parte sottostante per avere il plot delle sole linee
    %centrali
      
      %plot cono del sonar di prua (nero)
      c1t=plot3([pos(1) p1(i,1,1,1)],[pos(2) p1(i,2,1,1)],[pos(3) p1(i,3,1,1)],'k');
      c1p1=plot3(p1(i,1,1,1),p1(i,2,1,1),p1(i,3,1,1),'k*');
      c2t=plot3([pos(1) p1(i,1,1,2)],[pos(2) p1(i,2,1,2)],[pos(3) p1(i,3,1,2)],'k');
      c2p1=plot3(p1(i,1,1,2),p1(i,2,1,2),p1(i,3,1,2),'k*');
      c3t=plot3([pos(1) p1(i,1,1,3)],[pos(2) p1(i,2,1,3)],[pos(3) p1(i,3,1,3)],'k');
      c3p1=plot3(p1(i,1,1,3),p1(i,2,1,3),p1(i,3,1,3),'k*');
      c4t=plot3([pos(1) p1(i,1,2,1)],[pos(2) p1(i,2,2,1)],[pos(3) p1(i,3,2,1)],'k');
      c4p1=plot3(p1(i,1,2,1),p1(i,2,2,1),p1(i,3,2,1),'k*');
      c6t=plot3([pos(1) p1(i,1,2,3)],[pos(2) p1(i,2,2,3)],[pos(3) p1(i,3,2,3)],'k');
      c6p1=plot3(p1(i,1,2,3),p1(i,2,2,3),p1(i,3,2,3),'k*');
      c7t=plot3([pos(1) p1(i,1,3,1)],[pos(2) p1(i,2,3,1)],[pos(3) p1(i,3,3,1)],'k');
      c7p1=plot3(p1(i,1,3,1),p1(i,2,3,1),p1(i,3,3,1),'k*');
      c8t=plot3([pos(1) p1(i,1,3,2)],[pos(2) p1(i,2,3,2)],[pos(3) p1(i,3,3,2)],'k');
      c8p1=plot3(p1(i,1,3,2),p1(i,2,3,2),p1(i,3,3,2),'k*');
      c9t=plot3([pos(1) p1(i,1,3,3)],[pos(2) p1(i,2,3,3)],[pos(3) p1(i,3,3,3)],'k');
      c9p1=plot3(p1(i,1,3,3),p1(i,2,3,3),p1(i,3,3,3),'k*');
      
      %plot cono del sonar di destra (ciano)
      d1t=plot3([pos(1) p2(i,1,1,1)],[pos(2) p2(i,2,1,1)],[pos(3) p2(i,3,1,1)],'c');
      d1p1=plot3(p2(i,1,1,1),p2(i,2,1,1),p2(i,3,1,1),'c*');
      d2t=plot3([pos(1) p2(i,1,1,2)],[pos(2) p2(i,2,1,2)],[pos(3) p2(i,3,1,2)],'c');
      d2p1=plot3(p2(i,1,1,2),p2(i,2,1,2),p2(i,3,1,2),'c*');
      d3t=plot3([pos(1) p2(i,1,1,3)],[pos(2) p2(i,2,1,3)],[pos(3) p2(i,3,1,3)],'c');
      d3p1=plot3(p2(i,1,1,3),p2(i,2,1,3),p2(i,3,1,3),'c*');
      d4t=plot3([pos(1) p2(i,1,2,1)],[pos(2) p2(i,2,2,1)],[pos(3) p2(i,3,2,1)],'c');
      d4p1=plot3(p2(i,1,2,1),p2(i,2,2,1),p2(i,3,2,1),'c*');
      d6t=plot3([pos(1) p2(i,1,2,3)],[pos(2) p2(i,2,2,3)],[pos(3) p2(i,3,2,3)],'c');
      d6p1=plot3(p2(i,1,2,3),p2(i,2,2,3),p2(i,3,2,3),'c*');
      d7t=plot3([pos(1) p2(i,1,3,1)],[pos(2) p2(i,2,3,1)],[pos(3) p2(i,3,3,1)],'c');
      d7p1=plot3(p2(i,1,3,1),p2(i,2,3,1),p2(i,3,3,1),'c*');
      d8t=plot3([pos(1) p2(i,1,3,2)],[pos(2) p2(i,2,3,2)],[pos(3) p2(i,3,3,2)],'c');
      d8p1=plot3(p2(i,1,3,2),p2(i,2,3,2),p2(i,3,3,2),'c*');
      d9t=plot3([pos(1) p2(i,1,3,3)],[pos(2) p2(i,2,3,3)],[pos(3) p2(i,3,3,3)],'c');
      d9p1=plot3(p2(i,1,3,3),p2(i,2,3,3),p2(i,3,3,3),'c*');
      
      %plot cono del sonar di sinistra (verde)
      s1t=plot3([pos(1) p3(i,1,1,1)],[pos(2) p3(i,2,1,1)],[pos(3) p3(i,3,1,1)],'g');
      s1p1=plot3(p3(i,1,1,1),p3(i,2,1,1),p3(i,3,1,1),'g*');
      s2t=plot3([pos(1) p3(i,1,1,2)],[pos(2) p3(i,2,1,2)],[pos(3) p3(i,3,1,2)],'g');
      s2p1=plot3(p3(i,1,1,2),p3(i,2,1,2),p3(i,3,1,2),'g*');
      s3t=plot3([pos(1) p3(i,1,1,3)],[pos(2) p3(i,2,1,3)],[pos(3) p3(i,3,1,3)],'g');
      s3p1=plot3(p3(i,1,1,3),p3(i,2,1,3),p3(i,3,1,3),'g*');
      s4t=plot3([pos(1) p3(i,1,2,1)],[pos(2) p3(i,2,2,1)],[pos(3) p3(i,3,2,1)],'g');
      s4p1=plot3(p3(i,1,2,1),p3(i,2,2,1),p3(i,3,2,1),'g*');
      s6t=plot3([pos(1) p3(i,1,2,3)],[pos(2) p3(i,2,2,3)],[pos(3) p3(i,3,2,3)],'g');
      s6p1=plot3(p3(i,1,2,3),p3(i,2,2,3),p3(i,3,2,3),'g*');
      s7t=plot3([pos(1) p3(i,1,3,1)],[pos(2) p3(i,2,3,1)],[pos(3) p3(i,3,3,1)],'g');
      s7p1=plot3(p3(i,1,3,1),p3(i,2,3,1),p3(i,3,3,1),'g*');
      s8t=plot3([pos(1) p3(i,1,3,2)],[pos(2) p3(i,2,3,2)],[pos(3) p3(i,3,3,2)],'g');
      s8p1=plot3(p3(i,1,3,2),p3(i,2,3,2),p3(i,3,3,2),'g*');
      s9t=plot3([pos(1) p3(i,1,3,3)],[pos(2) p3(i,2,3,3)],[pos(3) p3(i,3,3,3)],'g');
      s9p1=plot3(p3(i,1,3,3),p3(i,2,3,3),p3(i,3,3,3),'g*');
    % Per avere il plot delle sole linee centrali, commentare fino a
    % questa riga
      
      pause(0.1) %permette di vedere il punto 'avanzare'
      
      %Cancella i punti relativi alle linee centrali ed al punto,
      %permettendone un avanzamento grafico
    delete(deep_purple); delete(c5t); delete(d5t); delete(s5t) ;delete(c5p1)
    delete(d5p); delete(s5p); delete(t1); delete(t2); delete(t3)
      
     % Cancella le rette dei coni dei sonar
      delete(c1t); delete(c2t); delete(c3t); delete(c4t);
      delete(c6t); delete(c7t); delete(c8t); delete(c9t);
      
      delete(d1t); delete(d2t); delete(d3t); delete(d4t);
      delete(d6t); delete(d7t); delete(d8t); delete(d9t);
      
      delete(s1t); delete(s2t); delete(s3t); delete(s4t);
      delete(s6t); delete(s7t); delete(s8t); delete(s9t);
      
      %cancella i punti di intersezioni delle rette dei coni
      delete(c1p1); delete(c2p1); delete(c3p1); delete(c4p1);
      delete(c6p1); delete(c7p1); delete(c8p1); delete(c9p1);
      
      delete(d1p1); delete(d2p1); delete(d3p1); delete(d4p1);
      delete(d6p1); delete(d7p1); delete(d8p1); delete(d9p1);
      
      delete(s1p1); delete(s2p1); delete(s3p1); delete(s4p1);
      delete(s6p1); delete(s7p1); delete(s8p1); delete(s9p1);
      
      hold on
end
