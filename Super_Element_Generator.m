function [Kc,Mc,Fc] = Super_Element_Generator(Nodes,IEN,normal_plane_vector,Geometry,Mat_Prop,Assigned_Type)

%% Model parameters
  Type = IEN(1:2,:)*0;
  Type(2,:) = 2;

  [ndof,nel,neq,nnp,LM] = Geometry_Model_Parameters_Preprocessor(Nodes,IEN,Type);


%% Extract material

  E = Mat_Prop.MoE(Assigned_Type)
  rho = Mat_Prop.rho(Assigned_Type)
  PoRat = Mat_Prop.PoRat(Assigned_Type)

  G = E/(2*(1+PoRat));

%% Section properties

  outer_radius=Geometry.Circular.Outer_Radius(Type(3,e));
  thickness=Geometry.Circular.Thickness(Type(3,e));

  r  = outer_radius;
  t  = thickness;
  ri = r - t;

  A  = pi*(r^2 - ri^2);

  Iy = pi*(r^4 - ri^4)/4;
  Iz = pi*(r^4 - ri^4)/4;
  Ix = pi*(r^4 - ri^4)/2;
  J  = pi*(r^4 - ri^4)/2;

  ky = 0.53;
  kz = 0.53;

%% Element matrices
  Ke_Total = zeros(12,12,nel);
  Me_Total = zeros(12,12,nel);

  for e = 1:nel

    X1 = Nodes(IEN(1,e),2:4);
    X2 = Nodes(IEN(2,e),2:4);
    [T,a] = Transformation_Matrix_Gen(X1,X2,normal_plane_vector);

    [KeT,MeT] = Truss_3D_Elemental_Stiffness_Matrix(E,rho,A,a);
    [KeB,MeB] = Beam_3D_Elemental_Stiffness_Matrix(G,E,rho,A,a,J,Ix,Iy,Iz,ky,kz);

    Ke_Total(:,:,e) = T'*(KeT+KeB)*T;
    Me_Total(:,:,e) = T'*(MeT+MeB)*T;

  end

%% Assembly
  K = Assembly_Truss(Ke_Total,Type,ndof,LM,nel,neq);
  M = Assembly_Truss(Me_Total,Type,ndof,LM,nel,neq);

  F = zeros(neq,1);

%% Static condensation
  bd = 1:12;
  id = 13:neq;

  Gg = -K(id,id)\K(id,bd);
  Tg = [eye(12); Gg];

  Kc = K(bd,bd) + K(bd,id)*Gg;
  Mc = Tg'*M*Tg;
  Fc = Gg'*F(id);

end
