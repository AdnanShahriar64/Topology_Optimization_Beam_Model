function [T,a] = Transformation_Matrix_Gen_Curved_Beam(Nodes,IEN,e,vref,ndofe)


  if nargin < 4
    vref = [0 0 1];
  end

  if nargin < 5
    ndofe = 12
  end



  X2=Nodes(IEN(2,e),2:4);
  X1=Nodes(IEN(1,e),2:4);

  d = (X2 - X1)';


  if e >1

    X2_P = Nodes(IEN(2,e-1),2:4);
    X1_P = Nodes(IEN(1,e-1),2:4);

    Start_vector = (X2_P - X1_P)';

    ez = cross(Start_vector/norm(Start_vector), d/norm(d))'

  end


  if e <size(IEN,2)

    X2_N = Nodes(IEN(2,e+1),2:4);
    X1_N = Nodes(IEN(1,e+1),2:4);

    End_vector = (X2_N - X1_N)';

    ez = cross(Start_vector/norm(Start_vector), d/norm(d))'

  end




  L = norm(d);


  a = L/2;

  ex = d/L;                                  % local x: along the beam

  ey = cross(ez,ex);                         % local y completes right-hand set

  R = [ex; ey; ez];                          % rows = local axes in global coords

  T=zeros(ndofe,ndofe);
  for i=1:4
      T(i*3-2:i*3,i*3-2:i*3)=R;
  end




end
