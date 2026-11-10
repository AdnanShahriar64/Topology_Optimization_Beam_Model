%function [T,a] = Transformation_Matrix_Gen_Curved_Beam(Nodes,IEN,e)

Nodes = [1  0     0     0;
         2  0     0.1  -1;
         3  0.5  -0.1  -1.1];


IEN =[1 2
      2 3];


  ndofe = 12

  X2=Nodes(IEN(2,e),2:4);
  X1=Nodes(IEN(1,e),2:4);

  d = (X2 - X1);
  L = norm(d);

  ey = d/L;                   % local y: along the beam


  Start_vector = [0 0 0];
  End_vector =  [0 0 0];
  ez_start = [0 0 0];
  ez_end = [0 0 0];


  if e >1

    X2_P = Nodes(IEN(2,e-1),2:4);
    X1_P = Nodes(IEN(1,e-1),2:4);

    Start_vector = (X2_P - X1_P);

    Start_vector = Start_vector/norm(Start_vector);

    ez_start = cross(Start_vector, ey);

  end

  if e <size(IEN,2)

    X2_N = Nodes(IEN(2,e+1),2:4);
    X1_N = Nodes(IEN(1,e+1),2:4);

    End_vector = (X2_N - X1_N);

    End_vector = End_vector/norm(End_vector);

    ez_end = cross(End_vector, ey);


  end


  ex = Start_vector;

  % ex = cross(ez,ey);                         % local y completes right-hand set

  R = transpose([ex; ey; ez]);                          % rows = local axes in global coords

  T=zeros(ndofe,ndofe);
  for i=1:4
      T(i*3-2:i*3,i*3-2:i*3)=R;
  end





