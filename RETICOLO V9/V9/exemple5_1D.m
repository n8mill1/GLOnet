%   1  D    exemple5_1D
%  VARIATION DE LA TRANSMISSION DANS l'ORDRE -1  FONCTION DE LAMBDA 
clear;
figure;
t=[];lambda=[];


teta0=0;nh=1;beta0=nh*sin(teta0*pi/180);
pol=-1; % 1:TE   -1:TM
parm=res0(pol);parm.not_io=1; 
parm.sym.x=0;% utilisation de la symetrie
nn=20;% ordres de fourier 

% description des textures 
textures{1}=1;   
textures{2}= 1.5; 
textures{3}={[-4,4],[1.5,1]  };


for LD=linspace(9,11,200);% longueur d'onde
D=10;% pas du reseau
teta0=0;nh=1;beta0=nh*sin(teta0*pi/180);
pol=-1; % 1:TE   -1:TM
parm=res0(pol);parm.not_io=1;   % initialisation des parametres par defaut
parm.sym.x=0;% utilisation de la symetrie
nn=10;% ordres de fourier 


% description des textures 
textures{1}=1;   
textures{2}= 1.5; 
textures{3}={[-4,4],[1.5,1]  };

% initialisation
aa=res1(LD,D,textures,nn,beta0,parm);
% definition du profil et calcul de la diffraction (on doit refaire ce calcul pour chaque LD)

profil={[0,20,0] ,[1,3,2]  };
    
ef=res2(aa,profil);

 % efficacite dans l'ordre -1 pour l'incident du haut 
t=[t,ef.inc_top_transmitted.efficiency{-1}];lambda=[lambda,LD];
plot(lambda,t);xlabel('longueur d''onde');ylabel('transmission dans l''ordre -1');title('VARIATION DE LA TRANSMISSION DANS l''ORDRE -1 FONCTION DE LAMBDA');pause(eps);
end





