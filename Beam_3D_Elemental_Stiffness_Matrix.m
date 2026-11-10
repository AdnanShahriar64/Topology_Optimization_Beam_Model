function[KeB,MeB]=Beam_3D_Elemental_Stiffness_Matrix(G,E,rho,A,a,J,Ix,Iy,Iz,ky,kz)


%% Obtain Stiffness matrix of Beam section for E=1, rho=1

% if e==1
%     fprintf('       Called Function (several times): Beam_3D_Elemental_Stiffness_Matrix\n')
%     fprintf('              Location:%s\n',pwd)
%     fprintf('              Purpose: Obtain Stiffness and Mass matrix for Bending component of each element\n\n')
% end

% | Input | Description | SI units |
% |---|---|---|
% | `G`   | Shear modulus, used for torsional stiffness and transverse shear stiffness. | Pa |
% | `Ey`  | Young’s modulus used for bending **about local \(y\)**, associated with displacement along local \(z\). | Pa |
% | `Ez`  | Young’s modulus used for bending **about local \(z\)**, associated with displacement along local \(y\). | Pa |
% | `rho` | Material mass density. | kg/m³ |
% | `A`   | Cross-sectional area. | m² |
% | `a`   | **Half the element length:** \(a=L/2\). | m |
% | `J`   | Saint-Venant torsion constant, used in torsional rigidity \(GJ\). Equals the polar area moment for circular sections. | m⁴ |
% | `Ix`  | Polar second moment of area about the longitudinal local \(x\)-axis, used for torsional rotary inertia: \(I_x=I_y+I_z\). | m⁴ |
% | `Iy`  | Second moment of area about local \(y\), used for bending rigidity \(E_yI_y\) and rotary inertia about \(y\). | m⁴ |
% | `Iz`  | Second moment of area about local \(z\), used for bending rigidity \(E_zI_z\) and rotary inertia about \(z\). | m⁴ |
% | `ky`  | Shear correction factor for transverse displacement along local \(y\); effective shear rigidity is \(k_yGA\). | Dimensionless |
% | `kz`  | Shear correction factor for transverse displacement along local \(z\); effective shear rigidity is \(k_zGA\). | Dimensionless |

Ey = E; Ez = E;
ax = a; ay = a; az = a;

V=A*a;
rx2=(Ix/A);

Phiz=3*Ez*Iz/(ky*G*A*az^2);
Phiy=3*Ey*Iy/(kz*G*A*ay^2);



%% Beam Stiffness part
    % TrK=zeros(12,12);
    %    u1     v1                              w1                              Ax1          Ay1                             Az1

   KeB1=[0      0                               0                               0            0                               0                               0    0                               0                               0           0                               0
         0      3*Ez*Iz/(2*az^3*(1+Phiz))       0                               0            0                               3*Ez*Iz/(2*az^2*(1+Phiz))       0   -3*Ez*Iz/(2*az^3*(1+Phiz))       0                               0           0                               3*Ez*Iz/(2*az^2*(1+Phiz))
         0      0                               3*Ey*Iy/(2*ay^3*(1+Phiy))       0           -3*Ey*Iy/(2*ay^2*(1+Phiy))       0                               0    0                              -3*Ey*Iy/(2*ay^3*(1+Phiy))       0          -3*Ey*Iy/(2*ay^2*(1+Phiy))       0
         0      0                               0                               G*J/(2*ax)   0                               0                               0    0                               0                              -G*J/(2*ax)  0                               0
         0      0                               0                               0            (4+Phiy)*Ey*Iy/(2*ay*(1+Phiy))  0                               0    0                               3*Ey*Iy/(2*ay^2*(1+Phiy))       0           (2-Phiy)*Ey*Iy/(2*ay*(1+Phiy))  0
         0      0                               0                               0            0                               (4+Phiz)*Ez*Iz/(2*az*(1+Phiz))  0   -3*Ez*Iz/(2*az^2*(1+Phiz))       0                               0           0                               (2-Phiz)*Ez*Iz/(2*az*(1+Phiz))
         0      0                               0                               0            0                               0                               0    0                               0                               0           0                               0
         0      0                               0                               0            0                               0                               0    3*Ez*Iz/(2*az^3*(1+Phiz))       0                               0           0                              -3*Ez*Iz/(2*az^2*(1+Phiz))
         0      0                               0                               0            0                               0                               0    0                               3*Ey*Iy/(2*ay^3*(1+Phiy))       0           3*Ey*Iy/(2*ay^2*(1+Phiy))       0
         0      0                               0                               0            0                               0                               0    0                               0                               G*J/(2*ax)  0                               0
         0      0                               0                               0            0                               0                               0    0                               0                               0           (4+Phiy)*Ey*Iy/(2*ay*(1+Phiy))  0
         0      0                               0                               0            0                               0                               0    0                               0                               0           0                               (4+Phiz)*Ez*Iz/(2*az*(1+Phiz))];


     KeB=KeB1+KeB1';
     for i=1:12
         KeB(i,i)=KeB1(i,i);
     end


%% Mass matrix
   Mz1=(78+147*Phiz+70*Phiz^2+63*Iz/(A*a^2))/(1+Phiz)^2;
   Mz2=((22+38.5*Phiz+17.5*Phiz^2)*a+(10.5-52.5*Phiz)*Iz/(A*a))/(1+Phiz)^2;
   Mz3=(27+63*Phiz+35*Phiz^2-63*Iz/(A*a^2))/(1+Phiz)^2;
   Mz4=(-(13+31.5*Phiz+17.5*Phiz^2)*a+(10.5-52.5*Phiz)*Iz/(A*a))/(1+Phiz)^2;
   Mz5=((8+14*Phiz+7*Phiz^2)*a^2+(28+35*Phiz+70*Phiz^2)*Iz/A)/(1+Phiz)^2;
   Mz6=(-(6+14*Phiz+7*Phiz^2)*a^2+(-7-35*Phiz+35*Phiz^2)*Iz/A)/(1+Phiz)^2;

   My1=(78+147*Phiy+70*Phiy^2+63*Iy/(A*a^2))/(1+Phiy)^2;
   My2=((22+38.5*Phiy+17.5*Phiy^2)*a+(10.5-52.5*Phiy)*Iy/(A*a))/(1+Phiy)^2;
   My3=(27+63*Phiy+35*Phiy^2-63*Iy/(A*a^2))/(1+Phiy)^2;
   My4=(-(13+31.5*Phiy+17.5*Phiy^2)*a+(10.5-52.5*Phiy)*Iy/(A*a))/(1+Phiy)^2;
   My5=((8+14*Phiy+7*Phiy^2)*a^2+(28+35*Phiy+70*Phiy^2)*Iy/A)/(1+Phiy)^2;
   My6=(-(6+14*Phiy+7*Phiy^2)*a^2+(-7-35*Phiy+35*Phiy^2)*Iy/A)/(1+Phiy)^2;

   MeB1=rho*V/105*[  0    0     0       0        0       0      0       0       0       0          0        0
                     0    Mz1   0       0        0       Mz2    0       Mz3     0       0          0        Mz4
                     0    0     My1     0       -My2     0      0       0       My3     0         -My4      0
                     0    0     0       70*rx2   0       0      0       0       0       35*rx2     0        0
                     0    0     0       0        My5     0      0       0       My4     0          My6      0
                     0    0     0       0        0       Mz5    0      -Mz4     0       0          0        Mz6
                     0    0     0       0        0       0      0       0       0       0          0        0
                     0    0     0       0        0       0      0       Mz1     0       0          0       -Mz2
                     0    0     0       0        0       0      0       0       My1     0          My2      0
                     0    0     0       0        0       0      0       0       0       70*rx2     0        0
                     0    0     0       0        0       0      0       0       0       0          My5      0
                     0    0     0       0        0       0      0       0       0       0          0        Mz5];


     MeB=MeB1+MeB1';
     for i=1:12
         MeB(i,i)=MeB1(i,i);
     end


end
