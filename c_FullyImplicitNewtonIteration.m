%%%%%%%%% inTheNameOfGod %%%%%%%%%
clc; clear;
tic
%% Gridding: Units in m %
n=41;
m=41;  
X=41*40/3.28; Y=41*40/3.28;
Z=1*25/3.28;
deltaX=40/3.28;
deltaZ=25/3.28;
deltaY=40/3.28;     %Grid numbers in y direction
% The gridding information is discussed in detail in the report. 
%% Defining rock and fluid properties
phi=0.15;
muo=1.2; muw=.8; %cp
cr=1e-6; %1/psi                    
co=20e-6;
cw=3e-6;
cf=co+cw;
Bob=1.5; Bwb=1; %STB/Bbl
P0=3600; %psia
initWaterSat=.15; %Not less than .2
resOilSat=.1;%
%% Defining computational specifics
ratio=0.02; %Differentiation Step
%% Reservoir schedule
bhp=1000; %psi
Qw=-800; %STBD
pi_1=36; pj_1=5;     %Producer #1 cordinates
pi_2=5; pj_2=5;      %Producer #2 cordinates
pi_3=36; pj_3=36;    %Producer #3 cordinates
pi_4=5; pj_4=36;     %Producer #4 cordinates
ii=21; ij=21;        %Injector cordinates
rw=0.1;
%% Converting units to SI 
muo=muo*1e-3; muw=muw*1e-3; %Pa-s
cr=cr/6894.76; cf=cf/6894.76; co=co/6894.76; cw=cw/6894.76; %1/pa
P0=P0*6894.76; bhp=bhp*6894.76; %Pa
P0INIT=P0*ones(m,n);
P0INITW=P0INIT-fitpcfunc(1.826e+06,-6.476,initWaterSat);
Qw=Qw*1.8403e-06; %m3/s 
qprimew=zeros(m,n);
qprimew(ii,ij)=Qw/deltaX/deltaY/deltaZ;
%% Defining K field%
load('K.mat');
%% Well treatment
WC=zeros(m,n);
re=.14*sqrt(deltaX^2+deltaY^2);
%Prod1
WC(pi_1,pj_1)=2*pi*K(pi_1,pj_1)*deltaZ/log(re/rw);
%Prod2
WC(pi_2,pj_2)=2*pi*K(pi_2,pj_2)*deltaZ/log(re/rw);
%Prod3
WC(pi_3,pj_3)=2*pi*K(pi_3,pj_3)*deltaZ/log(re/rw);
%Prod4
WC(pi_4,pj_4)=2*pi*K(pi_4,pj_4)*deltaZ/log(re/rw);
%Inj
WC(ii,ij)=2*pi*K(ii,ij)*deltaZ/log(re/rw);
%% Defining pressure, saturation, transmissibility/mobility, mobility, transmissibility matrices at time t
pressureO=P0*ones(m,n);
pressureW=pressureO-fitpcfunc(1.826e+06,-6.476,initWaterSat);
deltaPO=zeros(m,n);
deltaSW=zeros(m,n);
PO(2:m+1,2:n+1)=pressureO;
PW(2:m+1,2:n+1)=pressureW;
satw=initWaterSat*ones(m,n);
sato=(1-initWaterSat)*ones(m,n);
notTx=zeros(m,n+1);     %TransmissibilityX/mobility
notTy=zeros(m+1,n);     %TransmissibilityY/mobility
lambdaox=zeros(m,n+1);  lambdaoy=zeros(m+1,n);  % Oil transmissibilities in x and y direction
lambdawx=zeros(m,n+1);  lambdawy=zeros(m+1,n);  % Water transmissibilities in x and y direction
lambdao=zeros(m,n);
lambdaw=zeros(m,n);
Tox=zeros(m,n+1);   Toy=zeros(m+1,n);   % Oil mobilities in x and y direction
Twx=zeros(m,n+1);   Twy=zeros(m+1,n);   % Water mobilities in x and y direction 
for i=1:m
    for j=2:n
        notTx(i,j)=2/(deltaX*(deltaX/K(i,j-1)+deltaX/K(i,j)));
    end
end
for j=1:n
    for i=2:m
        notTy(i,j)=2/(deltaY*(deltaY/K(i-1,j)+deltaY/K(i,j)));
    end
end
%% Defining 3D matrices to store pressure and saturation at all the times
pressureO_3D=zeros(m,n,1);
pressureW_3D=zeros(m,n,1);
SATW=zeros(m,n,1);
SATO=zeros(m,n,1);
%% Initializing
pressureO_3D(:,:,1)=pressureO;   
pressureW_3D(:,:,1)=pressureW;
SATW(:,:,1)=satw;
SATO(:,:,1)=sato;
initOilInPlace=phi*X*Y*Z*(1-initWaterSat-resOilSat);
recoveryVector=zeros(1,1);
%% Solving algortithm
oilProdRate1Vector=zeros(1,1);
oilProdRate2Vector=zeros(1,1);
oilProdRate3Vector=zeros(1,1);
oilProdRate4Vector=zeros(1,1);
oilProdCum1Vector=zeros(1,1);
oilProdCum2Vector=zeros(1,1);
oilProdCum3Vector=zeros(1,1);
oilProdCum4Vector=zeros(1,1);
waterProdRate1Vector=zeros(1,1);
waterProdRate2Vector=zeros(1,1);
waterProdRate3Vector=zeros(1,1);
waterProdRate4Vector=zeros(1,1);
waterProdCum1Vector=zeros(1,1);
waterProdCum2Vector=zeros(1,1);
waterProdCum3Vector=zeros(1,1);
waterProdCum4Vector=zeros(1,1);
waterCut1Vector=zeros(1,1);
waterCut2Vector=zeros(1,1);
waterCut3Vector=zeros(1,1);
waterCut4Vector=zeros(1,1);
injWellBHPVector=zeros(1,1);
avOilPresVector=zeros(1,1);
avWaterSatVector=zeros(1,1);
timeVector=zeros(1,1);
deltaTVector=zeros(1,1);
A=zeros(2*m*n,2*m*n); B=zeros(2*m*n,1);
finalM=n+(m-1)*n;
t=1;
waterCut1=0; waterCut2=0; waterCut3=0; waterCut4=0;
TOL=600;
time=0;
Bo=Bob; Bw=Bwb;
%% Convergence Failure
error=false;
tsnumber=1;  
%%
while waterCut1<.95 || waterCut2<.95 || waterCut3<.95 || waterCut4<.95
    if error==false
        deltaT=tsnumber^3/2;
    else
        tsnumber=tsnumber-1;
        deltaT=tsnumber^3/2;
    end
    time=time+deltaT;
    fprintf('Time Step: %0f Seconds\n',deltaT)
    fprintf('Time: %0f Hours\n',time/3600)
    pressureOld=zeros(m,n);
    PO=zeros(m+2,n+2);
    PW=zeros(m+2,n+2);
    PO(2:m+1,2:n+1)=pressureO;
    PW(2:m+1,2:n+1)=pressureW;
    deltaPO=zeros(m,n);
    deltaSW=zeros(m,n);
    iter=1;
    sumNew=m*n*TOL+1;
    while sumNew>m*n*TOL
        sumOld=sumNew;
        fprintf('Iter No: %1g\n',iter)
        pressureO=PO(2:m+1,2:n+1);
        pressureW=PW(2:m+1,2:n+1);
        Cpoo=phi*(1-satw)*(cr+co)/Bo/deltaT;
        Cswo=-phi/Bo/deltaT;
        Cpow=phi*satw*(cr+cw)/Bw/deltaT;
        Csww=phi/Bw/deltaT-fitpcderfunc(1.826e+06,-6.476,satw).*Cpow;
        alpha=-Cswo./Csww;
        sWD=(satw-initWaterSat)/(1-initWaterSat-resOilSat);
        %% Calculating mobilities for internal grids 
        for i=2:m-1                         
            for j=2:n-1
                %% Oil mobilities
                if pressureO(i,j-1)>=pressureO(i,j)
                    lambdaox(i,j)=((1-sWD(i,j-1))^3.8)/muo/Bo;
                else
                    lambdaox(i,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                if pressureO(i+1,j)>=pressureO(i,j)
                    lambdaoy(i+1,j)=((1-sWD(i+1,j))^3.8)/muo/Bo;
                else
                    lambdaoy(i+1,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                %% Water Mobilities
                if pressureW(i,j-1)>=pressureW(i,j)
                    lambdawx(i,j)=(.4*sWD(i,j-1)^4)/muw/Bw;
                else
                    lambdawx(i,j)=(.4*sWD(i,j)^4)/muw/Bw;
                end
                if pressureW(i+1,j)>=pressureW(i,j)
                    lambdawy(i+1,j)=(.4*sWD(i+1,j)^4)/muw/Bw;
                else
                    lambdawy(i+1,j)=(.4*sWD(i,j)^4)/muw/Bw;
                end  
            end
        end
        %% Calculating mobilities for left face grids 
        for i=2:m-1                         
            for j=1
                %% Oil mobilities
                if pressureO(i+1,j)>=pressureO(i,j)
                    lambdaoy(i+1,j)=((1-sWD(i+1,j))^3.8)/muo/Bo;
                else
                    lambdaoy(i+1,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                %% Water mobilities
                if pressureW(i+1,j)>=pressureW(i,j)
                    lambdawy(i+1,j)=(.4*sWD(i+1,j)^4)/muw/Bw;
                else
                    lambdawy(i+1,j)=(.4*sWD(i,j)^4)/muw/Bw;
                end
            end
        end

        %% Calculating mobilities for right face grids 
        for i=2:m-1                         
            for j=n
                %% Oil mobilities
                if pressureO(i+1,j)>=pressureO(i,j)
                    lambdaoy(i+1,j)=((1-sWD(i+1,j))^3.8)/muo/Bo;
                else
                    lambdaoy(i+1,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                if pressureO(i,j)>=pressureO(i,j-1)
                    lambdaox(i,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                else
                    lambdaox(i,j)=((1-sWD(i,j-1))^3.8)/muo/Bo;
                end
                %% Water mobilities
                if pressureW(i+1,j)>=pressureW(i,j)
                    lambdawy(i+1,j)=(.4*sWD(i+1,j)^4)/muw/Bw;
                else
                    lambdawy(i+1,j)=(.4*sWD(i,j)^4)/muw/Bw;
                end
                if pressureW(i,j)>=pressureW(i,j-1)
                    lambdawx(i,j)=(.4*sWD(i,j)^4)/muw/Bw;
                else
                    lambdawx(i,j)=(.4*sWD(i,j-1)^4)/muw/Bw;
                end
            end
        end
        %% Calculating mobilities for top face grids 
        for i=1                         
            for j=2:n-1
                %% Oil mobilities
                if pressureO(i,j+1)>=pressureO(i,j)
                    lambdaox(i,j+1)=((1-sWD(i,j+1))^3.8)/muo/Bo;
                else
                    lambdaox(i,j+1)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                if pressureO(i+1,j)>=pressureO(i,j)
                    lambdaoy(i+1,j)=((1-sWD(i+1,j))^3.8)/muo/Bo;
                else
                    lambdaoy(i+1,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                %% Water mobilities
                if pressureW(i,j+1)>=pressureW(i,j)
                    lambdawx(i,j+1)=(.4*sWD(i,j+1)^4)/muw/Bw;
                else
                    lambdawx(i,j+1)=(.4*sWD(i,j)^4)/muw/Bw;
                end
                if pressureW(i+1,j)>=pressureW(i,j)
                    lambdawy(i+1,j)=(.4*sWD(i+1,j)^4)/muw/Bw;
                else
                    lambdawy(i+1,j)=(.4*sWD(i,j)^4)/muw/Bw;
                end
            end
        end
        %% Calculating mobilities for bottom face grids 
        for i=m                         
            for j=2:n-1
                %% Oil mobilities
                if pressureO(i,j-1)>=pressureO(i,j)
                    lambdaox(i,j)=((1-sWD(i,j-1))^3.8)/muo/Bo;
                else
                    lambdaox(i,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                end
                %% Water mobilities
                if pressureW(i,j-1)>=pressureW(i,j)
                    lambdawx(i,j)=(.4*sWD(i,j-1)^4)/muw/Bw;
                else
                    lambdawx(i,j)=(.4*sWD(i,j)^4)/muw/Bw;
                end
            end
        end
        %% Calculating mobilities for left top corner grid 
        i=1; j=1;
        % Oil mobilities
        if pressureO(i,j+1)>=pressureO(i,j)
            lambdaox(i,j+1)=((1-sWD(i,j+1))^3.8)/muo/Bo;
        else
            lambdaox(i,j+1)=((1-sWD(i,j))^3.8)/muo/Bo;
        end
        if pressureO(i+1,j)>=pressureO(i,j)
            lambdaoy(i+1,j)=((1-sWD(i+1,j))^3.8)/muo/Bo;
        else
            lambdaoy(i+1,j)=((1-sWD(i,j))^3.8)/muo/Bo;
        end
        % Water moblities
        if pressureW(i,j+1)>=pressureW(i,j)
            lambdawx(i,j+1)=(.4*sWD(i,j+1)^4)/muw/Bw;
        else
            lambdawx(i,j+1)=(.4*sWD(i,j)^4)/muw/Bw;
        end
        if pressureW(i+1,j)>=pressureW(i,j)
            lambdawy(i+1,j)=(.4*sWD(i+1,j)^4)/muw/Bw;
        else
            lambdawy(i+1,j)=(.4*sWD(i,j)^4)/muw/Bw;
        end
        %% Calculating mobilities for right top grid
        i=1; j=n;                       
        % Oil mobilities
        if pressureO(i+1,j)>=pressureO(i,j)
            lambdaoy(i+1,j)=((1-sWD(i+1,j))^3.8)/muo/Bo;
        else
            lambdaoy(i+1,j)=((1-sWD(i,j))^3.8)/muo/Bo;
        end
        % Water mobilities
        if pressureW(i+1,j)>=pressureW(i,j)
            lambdawy(i+1,j)=(.4*sWD(i+1,j)^4)/muw/Bw;
        else
            lambdawy(i+1,j)=(.4*sWD(i,j)^4)/muw/Bw;
        end
        Toy(i+1,j)=lambdaoy(i+1,j)*notTy(i+1,j);
        Twy(i+1,j)=lambdawy(i+1,j)*notTy(i+1,j);
        %% Calculating mobilities for bottom left corner grid
            % This has been calculated before
        %% Calculating mobilities for bottom right corner grid
        i=m; j=n;
        %	Oil mobilities
        if pressureO(i,j-1)>=pressureO(i,j)
            lambdaox(i,j)=((1-sWD(i,j-1))^3.8)/muo/Bo;
        else
            lambdaox(i,j)=((1-sWD(i,j))^3.8)/muo/Bo;
        end
        %	Water mobilities
        if pressureW(i,j-1)>=pressureW(i,j)
            lambdawx(i,j)=(.4*sWD(i,j-1)^4)/muw/Bw;
        else
            lambdawx(i,j)=(.4*sWD(i,j)^4)/muw/Bw;
        end
        %% Calculating transmissibilities in x direction
        for i=1:m
            for j=2:n
                Tox(i,j)=lambdaox(i,j)*notTx(i,j);
                Twx(i,j)=lambdawx(i,j)*notTx(i,j);
            end
        end
        %% Calculating transmissibilities in y direction
        for j=1:n
            for i=2:m
                Toy(i,j)=lambdaoy(i,j)*notTy(i,j);
                Twy(i,j)=lambdawy(i,j)*notTy(i,j);
            end
        end
        %% Calculating grid mobility
        for i=1:m
            for j=1:n
                lambdao(i,j)=((1-sWD(i,j))^3.8)/muo/Bo;
                lambdaw(i,j)=(.4*sWD(i,j)^4)/muw/Bw;
            end
        end
        %% Constructing A and B matrices for internal grids: A*X=B
        for i=2:m-1                         
            for j=2:n-1
                M=j+(i-1)*n;
                LOP=PO(i+1,j);    LWP=PW(i+1,j);
                TOP=PO(i,j+1);    TWP=PW(i,j+1);
                COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
                ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
                BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
                SW=satw(i,j);
                RondRwRondP_T=(RwCalc(LWP,(1+ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,(1-ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*TWP);
                RondRoRondP_T=(RoCalc(LOP,(1+ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,(1-ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*TOP);
                RondRwRondP_L=(RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*LWP);
                RondRoRondP_L=(RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*LOP);
                RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
                RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
                RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
                RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
                RondRwRondP_R=(RwCalc(LWP,TWP,CWP,(1+ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,(1-ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*RWP);
                RondRoRondP_R=(RoCalc(LOP,TOP,COP,(1+ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,(1-ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*ROP);
                RondRwRondP_B=(RwCalc(LWP,TWP,CWP,RWP,(1+ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,(1-ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*BWP);
                RondRoRondP_B=(RoCalc(LOP,TOP,COP,ROP,(1+ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,(1-ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*BOP);
                A(2*M-1,2*M-2*n-1)=RondRwRondP_T;
                A(2*M  ,2*M-2*n-1)=RondRoRondP_T;
                A(2*M-1,2*M-3)=RondRwRondP_L;
                A(2*M  ,2*M-3)=RondRoRondP_L;
                A(2*M-1,2*M-1)=RondRwRondP_C;
                A(2*M  ,2*M-1)=RondRoRondP_C;
                A(2*M-1,2*M)=RondRwRondSw_C;
                A(2*M  ,2*M)=RondRoRondSw_C;
                A(2*M-1,2*M+1)=RondRwRondP_R;
                A(2*M  ,2*M+1)=RondRoRondP_R;
                A(2*M-1,2*M+2*n-1)=RondRwRondP_B;
                A(2*M  ,2*M+2*n-1)=RondRoRondP_B;
                B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
                B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
            end
        end
        %% Constructing A and B matrices for left face grids: A*X=B
    for i=2:m-1                         
        for j=1
            M=j+(i-1)*n;
            LOP=PO(i+1,j);    LWP=PW(i+1,j);
            TOP=PO(i,j+1);    TWP=PW(i,j+1);
            COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
            ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
            BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
            SW=satw(i,j);
            RondRwRondP_T=(RwCalc(LWP,(1+ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,(1-ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*TWP);
            RondRoRondP_T=(RoCalc(LOP,(1+ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,(1-ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*TOP);
            RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
            RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
            RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRwRondP_R=(RwCalc(LWP,TWP,CWP,(1+ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,(1-ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*RWP);
            RondRoRondP_R=(RoCalc(LOP,TOP,COP,(1+ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,(1-ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*ROP);
            RondRwRondP_B=(RwCalc(LWP,TWP,CWP,RWP,(1+ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,(1-ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*BWP);
            RondRoRondP_B=(RoCalc(LOP,TOP,COP,ROP,(1+ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,(1-ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*BOP);
            A(2*M-1,2*M-2*n-1)=RondRwRondP_T;
            A(2*M  ,2*M-2*n-1)=RondRoRondP_T;
            A(2*M-1,2*M-1)=RondRwRondP_C;
            A(2*M  ,2*M-1)=RondRoRondP_C;
            A(2*M-1,2*M)=RondRwRondSw_C;
            A(2*M  ,2*M)=RondRoRondSw_C;
            A(2*M-1,2*M+1)=RondRwRondP_R;
            A(2*M  ,2*M+1)=RondRoRondP_R;
            A(2*M-1,2*M+2*n-1)=RondRwRondP_B;
            A(2*M  ,2*M+2*n-1)=RondRoRondP_B;
            B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
            B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        end
    end
    %% Constructing A and B matrices for right face grids: A*X=B
    for i=2:m-1                         
        for j=n
            M=j+(i-1)*n;
            LOP=PO(i+1,j);    LWP=PW(i+1,j);
            TOP=PO(i,j+1);    TWP=PW(i,j+1);
            COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
            ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
            BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
            SW=satw(i,j);
            RondRwRondP_T=(RwCalc(LWP,(1+ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,(1-ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*TWP);
            RondRoRondP_T=(RoCalc(LOP,(1+ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,(1-ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*TOP);
            RondRwRondP_L=(RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*LWP);
            RondRoRondP_L=(RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*LOP);
            RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
            RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
            RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRwRondP_B=(RwCalc(LWP,TWP,CWP,RWP,(1+ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,(1-ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*BWP);
            RondRoRondP_B=(RoCalc(LOP,TOP,COP,ROP,(1+ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,(1-ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*BOP);
            A(2*M-1,2*M-2*n-1)=RondRwRondP_T;
            A(2*M  ,2*M-2*n-1)=RondRoRondP_T;
            A(2*M-1,2*M-3)=RondRwRondP_L;
            A(2*M  ,2*M-3)=RondRoRondP_L;
            A(2*M-1,2*M-1)=RondRwRondP_C;
            A(2*M  ,2*M-1)=RondRoRondP_C;
            A(2*M-1,2*M)=RondRwRondSw_C;
            A(2*M  ,2*M)=RondRoRondSw_C;
            A(2*M-1,2*M+2*n-1)=RondRwRondP_B;
            A(2*M  ,2*M+2*n-1)=RondRoRondP_B;
            B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
            B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        end
    end
    %% Constructing A and B matrices for top face grids: A*X=B
    for i=1                         
        for j=2:n-1
            M=j+(i-1)*n;
            LOP=PO(i+1,j);    LWP=PW(i+1,j);
            TOP=PO(i,j+1);    TWP=PW(i,j+1);
            COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
            ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
            BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
            SW=satw(i,j);
            RondRwRondP_L=(RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*LWP);
            RondRoRondP_L=(RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*LOP);
            RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
            RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
            RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRwRondP_R=(RwCalc(LWP,TWP,CWP,(1+ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,(1-ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*RWP);
            RondRoRondP_R=(RoCalc(LOP,TOP,COP,(1+ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,(1-ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*ROP);
            RondRwRondP_B=(RwCalc(LWP,TWP,CWP,RWP,(1+ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,(1-ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*BWP);
            RondRoRondP_B=(RoCalc(LOP,TOP,COP,ROP,(1+ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,(1-ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*BOP);
            A(2*M-1,2*M-3)=RondRwRondP_L;
            A(2*M  ,2*M-3)=RondRoRondP_L;
            A(2*M-1,2*M-1)=RondRwRondP_C;
            A(2*M  ,2*M-1)=RondRoRondP_C;
            A(2*M-1,2*M)=RondRwRondSw_C;
            A(2*M  ,2*M)=RondRoRondSw_C;
            A(2*M-1,2*M+1)=RondRwRondP_R;
            A(2*M  ,2*M+1)=RondRoRondP_R;
            A(2*M-1,2*M+2*n-1)=RondRwRondP_B;
            A(2*M  ,2*M+2*n-1)=RondRoRondP_B;
            B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
            B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        end
    end
    %% Constructing A and B matrices for bottom face grids: A*X=B
    for i=m                         
        for j=2:n-1
            M=j+(i-1)*n;
            LOP=PO(i+1,j);    LWP=PW(i+1,j);
            TOP=PO(i,j+1);    TWP=PW(i,j+1);
            COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
            ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
            BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
            SW=satw(i,j);
            RondRwRondP_T=(RwCalc(LWP,(1+ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,(1-ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*TWP);
            RondRoRondP_T=(RoCalc(LOP,(1+ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,(1-ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*TOP);
            RondRwRondP_L=(RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*LWP);
            RondRoRondP_L=(RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*LOP);
            RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
            RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
            RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
            RondRwRondP_R=(RwCalc(LWP,TWP,CWP,(1+ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,(1-ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*RWP);
            RondRoRondP_R=(RoCalc(LOP,TOP,COP,(1+ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,(1-ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*ROP);
            A(2*M-1,2*M-2*n-1)=RondRwRondP_T;
            A(2*M  ,2*M-2*n-1)=RondRoRondP_T;
            A(2*M-1,2*M-3)=RondRwRondP_L;
            A(2*M  ,2*M-3)=RondRoRondP_L;
            A(2*M-1,2*M-1)=RondRwRondP_C;
            A(2*M  ,2*M-1)=RondRoRondP_C;
            A(2*M-1,2*M)=RondRwRondSw_C;
            A(2*M  ,2*M)=RondRoRondSw_C;
            A(2*M-1,2*M+1)=RondRwRondP_R;
            A(2*M  ,2*M+1)=RondRoRondP_R;
            B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
            B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        end
    end
        %% Constructing A and B matrices for left top corner grid: A*X=B
        i=1; j=1;
        M=j+(i-1)*n;
        LOP=PO(i+1,j);    LWP=PW(i+1,j);
        TOP=PO(i,j+1);    TWP=PW(i,j+1);
        COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
        ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
        BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
        SW=satw(i,j);
        RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
        RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
        RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRwRondP_R=(RwCalc(LWP,TWP,CWP,(1+ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,(1-ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*RWP);
        RondRoRondP_R=(RoCalc(LOP,TOP,COP,(1+ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,(1-ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*ROP);
        RondRwRondP_B=(RwCalc(LWP,TWP,CWP,RWP,(1+ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,(1-ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*BWP);
        RondRoRondP_B=(RoCalc(LOP,TOP,COP,ROP,(1+ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,(1-ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*BOP);
        A(2*M-1,2*M-1)=RondRwRondP_C;
        A(2*M  ,2*M-1)=RondRoRondP_C;
        A(2*M-1,2*M)=RondRwRondSw_C;
        A(2*M  ,2*M)=RondRoRondSw_C;
        A(2*M-1,2*M+1)=RondRwRondP_R;
        A(2*M  ,2*M+1)=RondRoRondP_R;
        A(2*M-1,2*M+2*n-1)=RondRwRondP_B;
        A(2*M  ,2*M+2*n-1)=RondRoRondP_B;
        B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
        B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        %% Constructing A and B matrices for right top corner grid: A*X=B
        i=1; j=n;                      
        M=j+(i-1)*n;
        LOP=PO(i+1,j);    LWP=PW(i+1,j);
        TOP=PO(i,j+1);    TWP=PW(i,j+1);
        COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
        ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
        BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
        SW=satw(i,j);
        RondRwRondP_L=(RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*LWP);
        RondRoRondP_L=(RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*LOP);
        RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
        RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
        RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRwRondP_B=(RwCalc(LWP,TWP,CWP,RWP,(1+ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,(1-ratio)*BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*BWP);
        RondRoRondP_B=(RoCalc(LOP,TOP,COP,ROP,(1+ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,(1-ratio)*BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*BOP);
        A(2*M-1,2*M-3)=RondRwRondP_L;
        A(2*M  ,2*M-3)=RondRoRondP_L;
        A(2*M-1,2*M-1)=RondRwRondP_C;
        A(2*M  ,2*M-1)=RondRoRondP_C;
        A(2*M-1,2*M)=RondRwRondSw_C;
        A(2*M  ,2*M)=RondRoRondSw_C;
        A(2*M-1,2*M+2*n-1)=RondRwRondP_B;
        A(2*M  ,2*M+2*n-1)=RondRoRondP_B;
        B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
        B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        %% Constructing A and B matrices for bottom left corner grid: A*X=B
        i=m; j=1;
        M=j+(i-1)*n;
        LOP=PO(i+1,j);    LWP=PW(i+1,j);
        TOP=PO(i,j+1);    TWP=PW(i,j+1);
        COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
        ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
        BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
        SW=satw(i,j);
        RondRwRondP_T=(RwCalc(LWP,(1+ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,(1-ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*TWP);
        RondRoRondP_T=(RoCalc(LOP,(1+ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,(1-ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*TOP);
        RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
        RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
        RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRwRondP_R=(RwCalc(LWP,TWP,CWP,(1+ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,(1-ratio)*RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*RWP);
        RondRoRondP_R=(RoCalc(LOP,TOP,COP,(1+ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,(1-ratio)*ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*ROP);
        A(2*M-1,2*M-2*n-1)=RondRwRondP_T;
        A(2*M  ,2*M-2*n-1)=RondRoRondP_T;
        A(2*M-1,2*M-1)=RondRwRondP_C;
        A(2*M  ,2*M-1)=RondRoRondP_C;
        A(2*M-1,2*M)=RondRwRondSw_C;
        A(2*M  ,2*M)=RondRoRondSw_C;
        A(2*M-1,2*M+1)=RondRwRondP_R;
        A(2*M  ,2*M+1)=RondRoRondP_R;
        B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
        B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        %% Constructing A and B matrices for bottom right corner grid: A*X=B
        i=m; j=n;
        M=j+(i-1)*n;
        LOP=PO(i+1,j);    LWP=PW(i+1,j);
        TOP=PO(i,j+1);    TWP=PW(i,j+1);
        COP=PO(i+1,j+1);  CWP=PW(i+1,j+1);
        ROP=PO(i+1,j+2);  RWP=PW(i+1,j+2);
        BOP=PO(i+2,j+1);  BWP=PW(i+2,j+1);
        SW=satw(i,j);
        RondRwRondP_T=(RwCalc(LWP,(1+ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,(1-ratio)*TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*TWP);
        RondRoRondP_T=(RoCalc(LOP,(1+ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,(1-ratio)*TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*TOP);
        RondRwRondP_L=(RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc((1+ratio)*LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*LWP);
        RondRoRondP_L=(RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc((1+ratio)*LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*LOP);
        RondRwRondP_C=(RwCalc(LWP,TWP,(1+ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,(1-ratio)*CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*CWP);
        RondRoRondP_C=(RoCalc(LOP,TOP,(1+ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,(1-ratio)*COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*COP);
        RondRwRondSw_C=(RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        RondRoRondSw_C=(RoCalc(LOP,TOP,COP,ROP,BOP,(1+ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)-RoCalc(LOP,TOP,COP,ROP,BOP,(1-ratio)*SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ))/(2*ratio*SW);
        A(2*M-1,2*M-2*n-1)=RondRwRondP_T;
        A(2*M  ,2*M-2*n-1)=RondRoRondP_T;
        A(2*M-1,2*M-3)=RondRwRondP_L;
        A(2*M  ,2*M-3)=RondRoRondP_L;
        A(2*M-1,2*M-1)=RondRwRondP_C;
        A(2*M  ,2*M-1)=RondRoRondP_C;
        A(2*M-1,2*M)=RondRwRondSw_C;
        A(2*M  ,2*M)=RondRoRondSw_C;
        B(2*M-1)=-RwCalc(LWP,TWP,CWP,RWP,BWP,COP,SW,i,j,t,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ);%-Rw;
        B(2*M)=-RoCalc(LOP,TOP,COP,ROP,BOP,SW,i,j,t,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ);%Ro;
        %% Calculating pressure 
        deltaPODeltaSWColumn=A\B;
        for u=1:m
            deltaPO(u,:)=deltaPODeltaSWColumn((u-1)*2*n+1:2:2*u*n-1,1);
        end
        PO(2:m+1,2:n+1)=PO(2:m+1,2:n+1)+deltaPO;
        for u=1:m
            deltaSW(u,:)=deltaPODeltaSWColumn((u-1)*2*n+2:2:2*u*n,1);
        end
        satw=satw+deltaSW;
        %%
        for i=2:m+1
            for j=2:n+1
                PW(i,j)=PO(i,j)-fitpcfunc(1.826e+06,-6.476,satw(i-1,j-1));
            end
        end
        sumNew=0;
        for i=1:m
            for j=1:n
                sumNew=sumNew+abs(deltaPO(i,j));
            end
        end
        if sumNew>sumOld && iter>1
            disp('Convergence Failed! Timestep Decreased!');
            time=time-deltaT;
            error=true;
            break
        else
            error=false;
        end
        pressureO=PO(2:m+1,2:n+1);
        pressureW=PW(2:m+1,2:n+1);
        iter=iter+1;
    end
    if error==true
        continue
    end
    error=false;
    fprintf('Total Number of Iterations in this Timestep: %1g\n',iter-1);
    %%
    pressureO_3D(:,:,t+1)=pressureO;
    for i=1:m
        for j=1:n
            pressureW(i,j)=pressureO(i,j)-fitpcfunc(1.826e+06,-6.476,satw(i,j));
        end
    end
    pressureW_3D(:,:,t+1)=pressureW(:,:);
    PressureW_3D_psi=pressureW_3D(:,:,t+1)/6894.76;
    pressureO_3D_psi=pressureO_3D(:,:,t+1)/6894.76;
    %% Calculating oil saturations
    for i=1:m
        for j=1:n
            sato(i,j)=1-satw(i,j);
        end
    end
    %% 
    SATW(:,:,t+1)=satw;
    SATO(:,:,t+1)=sato;
    qprimeo1=WC(pi_1,pj_1)*lambdao(pi_1,pj_1)*(pressureO(pi_1,pj_1)-bhp)/deltaX/deltaY/deltaZ;
    qprimeo2=WC(pi_2,pj_2)*lambdao(pi_2,pj_2)*(pressureO(pi_2,pj_2)-bhp)/deltaX/deltaY/deltaZ;
    qprimeo3=WC(pi_3,pj_3)*lambdao(pi_3,pj_3)*(pressureO(pi_3,pj_3)-bhp)/deltaX/deltaY/deltaZ;
    qprimeo4=WC(pi_4,pj_4)*lambdao(pi_4,pj_4)*(pressureO(pi_3,pj_3)-bhp)/deltaX/deltaY/deltaZ;
    Qprimeo1=qprimeo1*deltaX*deltaY*deltaZ;
    Qprimeo2=qprimeo2*deltaX*deltaY*deltaZ;
    Qprimeo3=qprimeo3*deltaX*deltaY*deltaZ;
    Qprimeo4=qprimeo4*deltaX*deltaY*deltaZ;
    Qprimew1=WC(pi_1,pj_1)*(lambdaw(pi_1,pj_1))*(pressureW(pi_1,pj_1)-bhp);
    Qprimew2=WC(pi_2,pj_2)*(lambdaw(pi_2,pj_2))*(pressureW(pi_2,pj_2)-bhp);
    Qprimew3=WC(pi_3,pj_3)*(lambdaw(pi_3,pj_3))*(pressureW(pi_3,pj_3)-bhp);
    Qprimew4=WC(pi_4,pj_4)*(lambdaw(pi_4,pj_4))*(pressureW(pi_4,pj_4)-bhp);
    waterCut1=Qprimew1/(Qprimew1+Qprimeo1);
    waterCut2=Qprimew2/(Qprimew2+Qprimeo2);
    waterCut3=Qprimew3/(Qprimew3+Qprimeo3);
    waterCut4=Qprimew4/(Qprimew4+Qprimeo4);
    fprintf('Water Cut #1: %0f\n',waterCut1)
    fprintf('Water Cut #2: %0f\n',waterCut2)
    fprintf('Water Cut #3: %0f\n',waterCut3)
    fprintf('Water Cut #4: %0f\n',waterCut4)
    timeVector(t+1)=time;
    deltaTVector(t+1)=deltaT;
    oilProdRate1Vector(t+1)=Qprimeo1;
    oilProdRate2Vector(t+1)=Qprimeo2;
    oilProdRate3Vector(t+1)=Qprimeo3;
    oilProdRate4Vector(t+1)=Qprimeo4;
    oilProdCum1Vector(t+1)=oilProdCum1Vector(t)+((oilProdRate1Vector(t)+oilProdRate1Vector(t+1))*deltaT)/2;
    oilProdCum2Vector(t+1)=oilProdCum2Vector(t)+((oilProdRate2Vector(t)+oilProdRate2Vector(t+1))*deltaT)/2;
    oilProdCum3Vector(t+1)=oilProdCum3Vector(t)+((oilProdRate3Vector(t)+oilProdRate3Vector(t+1))*deltaT)/2;
    oilProdCum4Vector(t+1)=oilProdCum4Vector(t)+((oilProdRate4Vector(t)+oilProdRate4Vector(t+1))*deltaT)/2;
    waterProdRate1Vector(t+1)=Qprimew1;
    waterProdRate2Vector(t+1)=Qprimew2;
    waterProdRate3Vector(t+1)=Qprimew3;
    waterProdRate4Vector(t+1)=Qprimew4;
    waterProdCum1Vector(t+1)=waterProdCum1Vector(t)+((waterProdRate1Vector(t)+waterProdRate1Vector(t+1))*deltaT)/2;
    waterProdCum2Vector(t+1)=waterProdCum2Vector(t)+((waterProdRate2Vector(t)+waterProdRate2Vector(t+1))*deltaT)/2;
    waterProdCum3Vector(t+1)=waterProdCum3Vector(t)+((waterProdRate3Vector(t)+waterProdRate3Vector(t+1))*deltaT)/2;
    waterProdCum4Vector(t+1)=waterProdCum4Vector(t)+((waterProdRate4Vector(t)+waterProdRate4Vector(t+1))*deltaT)/2;
    recoveryVector(t+1)=(oilProdCum1Vector(t+1)+oilProdCum2Vector(t+1)+oilProdCum3Vector(t+1)+oilProdCum4Vector(t+1))/initOilInPlace*Bo;
    lambdaoInj=((1-sWD(ii,ij))^3.8)/muo/Bo;
    lambdawInj=(.4*sWD(ii,ij)^4)/muw/Bw;
    injWellBHPVector(t+1)=pressureW(ii,ij)-Qw/WC(ii,ij)/(Bo/Bw*lambdaoInj+lambdawInj);
    waterCut1Vector(t+1)=waterCut1;
    waterCut2Vector(t+1)=waterCut2;
    waterCut3Vector(t+1)=waterCut3;
    waterCut4Vector(t+1)=waterCut4;
    dintpdxdy=0;
    dintswdxdy=0;
    for i=1:m
        for j=1:n
            dintpdxdy=dintpdxdy+pressureO(i,j)*deltaX*deltaY;
            dintswdxdy=dintswdxdy+satw(i,j)*deltaX*deltaY;
        end
    end
    avOilPres=dintpdxdy/X/Y;
    avWaterSat=dintswdxdy/X/Y;
    avOilPresVector(t+1)=avOilPres;
    avWaterSatVector(t+1)=avWaterSat;
    fprintf('Recovery: %0f\n\n',recoveryVector(t+1))
    t=t+1;
    tsnumber=tsnumber+1;
    figure(1)
    movegui(figure(1),'west');
    imagesc(satw);
    axis equal
    axis tight
    colorbar
    str1 = sprintf('Water Saturation at T=%f Hours', time/3600);
    title(str1);
    figure(2)
    imagesc(pressureO/6894.76),
    axis equal
    axis tight
    colorbar
    str2 = sprintf('Oil Pressure (PSI) at T=%f Hours', time/3600);
    title(str2);
    drawnow
end
%% Almost there!
recoveryVectorSize=size(recoveryVector);
ultimateRecoveryFactor=recoveryVector(recoveryVectorSize(2));
timeVector=timeVector/3600;
oilProdRate1Vector=oilProdRate1Vector*543439.65;
oilProdRate2Vector=oilProdRate2Vector*543439.65;
oilProdRate3Vector=oilProdRate3Vector*543439.65;
oilProdRate4Vector=oilProdRate4Vector*543439.65;
waterProdRate1Vector=waterProdRate1Vector*543439.65;
waterProdRate2Vector=waterProdRate2Vector*543439.65;
waterProdRate3Vector=waterProdRate3Vector*543439.65;
waterProdRate4Vector=waterProdRate4Vector*543439.65;
oilProdCum1Vector=oilProdCum1Vector*6.28981;
oilProdCum2Vector=oilProdCum2Vector*6.28981;
oilProdCum3Vector=oilProdCum3Vector*6.28981;
oilProdCum4Vector=oilProdCum4Vector*6.28981;
waterProdCum1Vector=waterProdCum1Vector*6.28981;
waterProdCum2Vector=waterProdCum2Vector*6.28981;
waterProdCum3Vector=waterProdCum3Vector*6.28981;
waterProdCum4Vector=waterProdCum4Vector*6.28981;
injWellBHPVector=injWellBHPVector/6894.76;
avOilPresVector=avOilPresVector/6894.76;
%% The Results :)
fprintf('Ultimate Recovery Factor: %0f\n\n',ultimateRecoveryFactor)
figure(3)
plot(timeVector,oilProdRate1Vector,timeVector,oilProdRate2Vector,timeVector,oilProdRate3Vector,timeVector,oilProdRate4Vector,timeVector,waterProdRate1Vector,timeVector,waterProdRate2Vector,timeVector,waterProdRate3Vector,timeVector,waterProdRate4Vector);
title('Oil and Water Production Rates');
ylabel('Water and Oil Production Rates (STB/Day)'); xlabel('Time (Hours)');
legend('Oil #1','Oil #2','Oil #3','Oil #4','Water #1','Water #2','Water #3','Water #4')
figure(4)
plot(timeVector,oilProdCum1Vector,timeVector,oilProdCum2Vector,timeVector,oilProdCum3Vector,timeVector,oilProdCum4Vector,timeVector,waterProdCum1Vector,timeVector,waterProdCum2Vector,timeVector,waterProdCum3Vector,timeVector,waterProdCum4Vector);
title('Oil and Water Cumulative Production');
ylabel('Water and Oil Cumulative Production (STB)'); xlabel('Time (Hours)');
legend('Oil #1','Oil #2','Oil #3','Oil #4','Water #1','Water #2','Water #3','Water #4')
figure(5)
plot(timeVector,injWellBHPVector);
title('Injection Well Bottomhole Pressure');
ylabel('Pressure (PSI)'); xlabel('Time (Hours)');
figure(6)
plot(timeVector,recoveryVector);
title('Recovery Factor');
ylabel('Recovery Factor'); xlabel('Time (Hours)');
figure(7)
plot(timeVector,avOilPresVector);
title('Average Oil Pressure');
ylabel('Average Oil Pressure (PSI)'); xlabel('Time (Hours)');
figure(8)
plot(timeVector,avWaterSatVector);
title('Average Water Saturation');
ylabel('Average Water Saturation'); xlabel('Time (Hours)');
figure(9)
plot(timeVector,waterCut1Vector,timeVector,waterCut2Vector,timeVector,waterCut3Vector,timeVector,waterCut4Vector);
title('Water Cut');
legend('#1','#2','#3','#4')
figure(10)
plot(timeVector,deltaTVector);
title('Adaptive Timestep (s) vs. Try No.');
toc

%% Functions
function  Ro = RoCalc (LOP,TOP,COP,ROP,BOP,SW,ii,jj,tt,SATW,pressureO_3D,Cswo,Cpoo,Tox,Toy,bhp,WC,lambdao,deltaX,deltaY,deltaZ)
   Ro=Cpoo(ii,jj)*(COP-pressureO_3D(ii,jj,tt))+Cswo*(SW-SATW(ii,jj,tt))...
        -Tox(ii,jj+1)*(ROP-COP)+Tox(ii,jj)*(COP-LOP)...
        -Toy(ii,jj)*(TOP-COP)+Toy(ii+1,jj)*(COP-BOP)+WC(ii,jj)*lambdao(ii,jj)*(COP-bhp)/deltaX/deltaY/deltaZ;
end

function  Rw = RwCalc (LWP,TWP,CWP,RWP,BWP,COP,SW,ii,jj,tt,SATW,pressureO_3D,Csww,Cpow,Twx,Twy,bhp,WC,qprimew,lambdaw,deltaX,deltaY,deltaZ)
   Rw=Cpow(ii,jj)*(COP-pressureO_3D(ii,jj,tt))+Csww(ii,jj)*(SW-SATW(ii,jj,tt))...
        -Twx(ii,jj+1)*(RWP-CWP)+Twx(ii,jj)*(CWP-LWP)...
        -Twy(ii,jj)*(TWP-CWP)+Twy(ii+1,jj)*(CWP-BWP)+qprimew(ii,jj)+WC(ii,jj)*lambdaw(ii,jj)*(CWP-bhp)/deltaX/deltaY/deltaZ;
end

function son = fitpcfunc(a,b,x)
    son=a*exp(b*x);
end

function daughter = fitpcderfunc(a,b,x)
    daughter=a*b*exp(b*x);
end