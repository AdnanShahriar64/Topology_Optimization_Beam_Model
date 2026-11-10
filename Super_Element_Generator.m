function [K] = Super_Element_Generator(Nodes,IEN,E,rho,PoRat,outer_radius,thickness)

% IEN= IEN_L;
% Nodes=Nodes_L;
% E = Mat_Prop.MoE;
% rho = Mat_Prop.rho;
% PoRat = Mat_Prop.PoRat;
% outer_radius = 0.001;
% thickness = 0.0001;

  Type = IEN(1:2,:)*0;
  Type(2,:)=2;

  [ndof,nel,neq,nnp,LM]=Geometry_Model_Parameters_Preprocessor(Nodes,IEN,Type);

  G=E/(2*(1+Mat_Prop.PoRat));

% Get geometry data

  r = outer_radius;
  t = thickness;

  A=pi*(r^2-(r-t)^2);

  Iz=(pi*r^4)/4;
  Iy=(pi*r^4)/4;
  Ix=(pi*r^4)/4;

  J=pi*(2*r)^4/32;

  ky = 10^6;
  kz = 10^6;


%% Truss section
  KeT=zeros(12,12,nel);
  MeT=KeT;

  %% Beam section
  KeB=zeros(12,12,nel);
  MeB=KeB;

  for e=1:nel

    X1 = Nodes(IEN(1,e),2:4);
    X2 = Nodes(IEN(2,e),2:4);
    [T,a] = Transformation_Matrix_Gen(X1,X2);

    % Get local elemental stiffness matrix for truss with E_Rat=1, rho=1
      [KeT,MeT]=Truss_3D_Elemental_Stiffness_Matrix(E,rho,A,a);

    % Get local elemental stiffness matrix for Beam with E_Rat=1, rho=1
      [KeB,MeB]=Beam_3D_Elemental_Stiffness_Matrix(G,E,rho,A,a,J,Ix,Iy,Iz,ky,kz);



       Ke = KeT+KeB;
       Me = MeT+MeB;

       Ke_Total(:,:,e) = T'*Ke*T;
       Me_Total(:,:,e) = T'*Me*T;


  end


  [K] = Assembly_Truss(Ke_Total,Type,ndof,LM,nel,neq);
  [M] = Assembly_Truss(Me_Total,Type,ndof,LM,nel,neq);

  F = zeros(neq,1);

  bd = 1:12;
  id = 13:neq;

  Gg  = -K(id,id)\K(id,bd);
  Tg = [eye(12); Gg];

  Kc = K(bd,bd) + K(bd,id)*Gg;
  Mc = Tg'*M*Tg;
  % Fc = F(bd) + Gg'*F(id);
  Fc = Gg'*F(id);


 end







