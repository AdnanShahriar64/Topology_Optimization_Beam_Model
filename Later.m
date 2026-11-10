function [K] = Super_Element_Generator(Nodes,IEN,Type,Mat_Prop,Geometry)


  [ndof,nel,neq,nnp,LM]=Geometry_Model_Parameters_Preprocessor(Nodes,IEN,Type);


  E = Mat_Prop.MoE(Type(4,e))
  rho = Mat_Prop.rho(Type(4,e));
  G=E/(2*(1+Mat_Prop.PoRat(Type(4,e))));

% Get geometry data



  Iz=(pi*r^4)/4;
  Iy=(pi*r^4)/4;
  Ix=(pi*r^4)/4;

  J=pi*(2*r)^4/32;


%% Truss section
  KeT=zeros(12,12,nel);
  MeT=KeT;

  %% Beam section
  KeB=zeros(12,12,nel);
  MeB=KeB;

  for e=1:nel


    r = Geometry.Circular.Outer_Radius(Type(3,e);
    t = Geometry.Circular.Thickness(Type(3,e);

    A=pi*(r^2-(r-t)^2);

    X1 = Nodes(IEN(1,e),2:4);
    X2 = Nodes(IEN(2,e),2:4);
    [T,a] = Transformation_Matrix_Gen(X1,X2);

    % Get local elemental stiffness matrix for truss with E_Rat=1, rho=1
      [KeT,MeT]=Truss_3D_Elemental_Stiffness_Matrix(E,rho,A,a);

    % Get local elemental stiffness matrix for Beam with E_Rat=1, rho=1
      [KeB,MeB]=Beam_3D_Elemental_Stiffness_Matrix(G,E,rho,A,a,J,Ix,Iy,Iz,ky,kz)



         Ke = KeT+KeB;

         Ke_Total(:,:,e) = T'*Ke*T;


  end


  [K] = Assembly_Truss(Ke_Total,Type,ndof,LM,E_List,neq);





  end







