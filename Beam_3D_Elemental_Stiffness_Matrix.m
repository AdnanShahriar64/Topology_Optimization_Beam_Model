function[KeB,MeB]=Beam_3D_Elemental_Stiffness_Matrix(A,a,ax,ay,az,G,Ey,Ez,J,Ix,Iy,Iz)


%% Obtain Stiffness matrix of Beam section for E=1, rho=1

% if e==1
%     fprintf('       Called Function (several times): Beam_3D_Elemental_Stiffness_Matrix\n')
%     fprintf('              Location:%s\n',pwd)
%     fprintf('              Purpose: Obtain Stiffness and Mass matrix for Truss section for each element\n\n') 
% end


E=1;
rho=1;


V=A*a;
rx2=(Ix/A);

%% Truss part
%     AxK=zeros(12,12);
%     AxK(1,1)=A*E/(2*a);
%     AxK(1,6+1)=-A*E/(2*a);
%     AxK(6+1,1)=-A*E/(2*a);
%     AxK(6+1,6+1)=A*E/(2*a);
    
%% Beam Stiffness part
    % TrK=zeros(12,12);
    %    u1     v1                  w1                  Ax1          Ay1                 Az1
               
   KeB1=[0      0                   0                   0            0                   0                   0    0                   0                  0           0                   0
         0      3*Ez*Iz/(2*az^3)    0                   0            0                   3*Ez*Iz/(2*az^2)    0   -3*Ez*Iz/(2*az^3)    0                  0           0                   3*Ez*Iz/(2*az^2)
         0      0                   3*Ey*Iy/(2*ay^3)    0           -3*Ey*Iy/(2*ay^2)    0                   0    0                  -3*Ey*Iy/(2*ay^3)   0          -3*Ey*Iy/(2*ay^2)    0
         0      0                   0                   G*J/(2*ax)   0                   0                   0    0                   0                 -G*J/(2*ax)   0                   0
         0      0                   0                   0            2*Ey*Iy/ay          0                   0    0                   3*Ey*Iy/(2*ay^2)   0           Ey*Iy/ay            0                    
         0      0                   0                   0            0                   2*Ez*Iz/az          0   -3*Ez*Iz/(2*az^2)    0                  0           0                   Ez*Iz/az               
         0      0                   0                   0            0                   0                   0    0                   0                  0           0                   0                      
         0      0                   0                   0            0                   0                   0    3*Ez*Iz/(2*az^3)    0                  0           0                  -3*Ez*Iz/(2*az^2)
         0      0                   0                   0            0                   0                   0    0                   3*Ey*Iy/(2*ay^3)   0           3*Ey*Iy/(2*ay^2)    0                      
         0      0                   0                   0            0                   0                   0    0                   0                  G*J/(2*ax)  0                   0                      
         0      0                   0                   0            0                   0                   0    0                   0                  0           2*Ey*Iy/ay          0
         0      0                   0                   0            0                   0                   0    0                   0                  0           0                   2*Ez*Iz/az         ];

    
     KeB=KeB1+KeB1';
     for i=1:12
         KeB(i,i)=KeB1(i,i);
     end
     
     KeB=E*KeB;
    
     
     % KeI=TrK+AxK;
     
%% Mass matrix
   MeB1=rho*V/105*[ 0    0     0       0        0       0      0       0       0       0          0        0
                     0    78    0       0        0       22*a   0       27      0       0          0       -13*a
                     0    0     78      0       -22*a    0      0       0       27      0          13*a     0
                     0    0     0       70*rx2   0       0      0       0       0       35*rx2     0        0
                     0    0     0       0        8*a^2   0      0       0      -13*a    0         -6*a^2    0
                     0    0     0       0        0       8*a^2  0       13*a    0       0          0       -6*a^2
                     0    0     0       0        0       0      0       0       0       0          0        0   
                     0    0     0       0        0       0      0       78      0       0          0       -22*a
                     0    0     0       0        0       0      0       0       78      0          22*a     0
                     0    0     0       0        0       0      0       0       0       70*rx2     0        0
                     0    0     0       0        0       0      0       0       0       0          8*a^2    0
                     0    0     0       0        0       0      0       0       0       0          0        8*a^2];
   
                 
     MeB=MeB1+MeB1';
     for i=1:12
         MeB(i,i)=MeB1(i,i);
     end      
   
%      MeAx=zeros(12,12);
%      MeAx(1,1)=70;
%      MeAx(1,6+1)=35;
%      MeAx(6+1,1)=35;
%      MeAx(6+1,6+1)=70;
%      MeAx=MeAx*(rho*V)/105;
     
     %MeI=MeAx+MeB;
     
end
