function[KeT,MeT]=Truss_3D_Elemental_Stiffness_Matrix(A,a)


%% Obtain Stiffness matrix of Truss section for E=1, rho=1

% if e==1
%     fprintf('       Called Function (several times): Truss_3D_Elemental_Stiffness_Matrix\n')
%     fprintf('              Location:%s\n',pwd)
%     fprintf('              Purpose: Obtain Stiffness and Mass matrix for Truss section for each element\n') 
% end
    V=A*a;

    E=1;
    rho=1;
    
%% Truss part
    KeT=zeros(12,12);
    KeT(1,1)=A*E/(2*a);
    KeT(1,6+1)=-A*E/(2*a);
    KeT(6+1,1)=-A*E/(2*a);
    KeT(6+1,6+1)=A*E/(2*a);
    
     MeT=zeros(12,12);
     MeT(1,1)=70;
     MeT(1,6+1)=35;
     MeT(6+1,1)=35;
     MeT(6+1,6+1)=70;
     MeT=MeT*(rho*V)/105;
    




end