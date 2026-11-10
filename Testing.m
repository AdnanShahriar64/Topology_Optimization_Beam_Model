
Nodes = [1  0     0     0;
         2  0     0.1  -1;
         3  0.5  -0.1  -1.1
         4  0.5  -0.1   0.1];


IEN =[1 2 3
      2 3 4];


Type = [1 2 1;  % 1 standard beam, 2 curved beam
        2 2 2;  % No of nodes
        1 1 1;  % Group
        1 1 2]; % Material type, may be useful for joint approximation

nel = size(IEN,2);

Geometry.Circular.Outer_Radius =     0.023622;
Geometry.Circular.Thickness = 0.0012446;


Mat_Prop.MoE = 10^9;
Mat_Prop.PoRat = 0.3;
Mat_Prop.rho = 1000;


NB      = zeros(4,nel);
for e = 1:nel

  if Type(1,e) == 2

    n1 = IEN(1,e);
    n2 = IEN(2,e);

    n1_far = 0;
    n2_far = 0;
    n1_count = 0;
    n2_count = 0;
    for i = 1:nel
        if i ~= e
            if IEN(1,i) == n1
                n1_far = IEN(2,i);
                n1_count = n1_count + 1;
            elseif IEN(2,i) == n1
                n1_far = IEN(1,i);
                n1_count = n1_count + 1;
            end

            if IEN(1,i) == n2
                n2_far = IEN(2,i);
                n2_count = n2_count + 1;
            elseif IEN(2,i) == n2
                n2_far = IEN(1,i);
                n2_count = n2_count + 1;
            end
        end
    end

    NB(1,e) = n1_far;
    NB(2,e) = n1;
    NB(3,e) = n2;
    NB(4,e) = n2_far;

  end

end





  e = 2;

  ndofe = 12

  X0 = Nodes(NB(1,e),2:4);
  X1 = Nodes(NB(2,e),2:4);
  X2 = Nodes(NB(3,e),2:4);
  X3 = Nodes(NB(4,e),2:4);

  c  = X2 - X1;
  Lc = norm(c);
  ch = c/Lc;

  t1 = (X1 - X0)/norm(X1 - X0);
  t2 = (X3 - X2)/norm(X3 - X2);

  t2c = 2*dot(t1,ch)*ch - t1;
  err = norm(t2 - t2c);




  C = (X1 + X2)/2;
  R = Lc/2;

  b1 = cross(t1,ch);
  b1 = b1/norm(b1);

  b2 = cross(ch,t2);
  b2 = b2/norm(b2);

  b =  b1 + b2;

  lead = 0;

  if norm(b1-b2)>10^-5

    fprintf("2 connecting elements are not parallel with angle =%f degrees\n", acos(dot(b1,b2))*180/pi);

    if lead == 1
      b = b1;
      fprintf("Set previous element as leading element\n");
    elseif lead ==2
      b = b2;
      fprintf("Set Next element as leading element\n");
    end

  end







  if norm(b) < 1e-8
      error('Tangents along the chord: plane undefined.');
  end
  b = b/norm(b);

  t1c = cross(ch,b);
  ang1 = acos(min(1,dot(t1, t1c)));
  ang2 = acos(min(1,dot(t2,-t1c)));




nArc = 16;

b1 = cross(t1,ch);
b2 = cross(ch,t2);

Nodes_L = zeros(nArc+1,4);
Nodes_L(1,:) = [1 X1];
Nodes_L(2,:) = [2 X2];
for i = 1:nArc-1
    phi = pi*i/nArc;
    Nodes_L(i+2,1)   = i+2;
    Nodes_L(i+2,2:4) = C - R*cos(phi)*ch + R*sin(phi)*t1c;
end

IEN_L = zeros(2,nArc);
IEN_L(1,1) = 1;
for i = 1:nArc-1
    IEN_L(2,i)   = i+2;
    IEN_L(1,i+1) = i+2;
end
IEN_L(2,nArc) = 2;

outer_radius=Geometry.Circular.Outer_Radius(Type(3,e));
thickness=Geometry.Circular.Thickness(Type(3,e));

[Ke,Me] = Super_Element_Generator(Nodes_L,IEN_L,Mat_Prop.MoE,Mat_Prop.rho,Mat_Prop.PoRat,outer_radius,thickness);







figure
ord = [1 3:nArc+1 2];
plot3(Nodes_L(ord,2),Nodes_L(ord,3),Nodes_L(ord,4))
hold on
scatter3(Nodes_L(:,2),Nodes_L(:,3),Nodes_L(:,4))
plot3([X0(1) X1(1)],[X0(2) X1(2)],[X0(3) X1(3)],'k')
plot3([X2(1) X3(1)],[X2(2) X3(2)],[X2(3) X3(3)],'k')

quiver3(X1(1),X1(2),X1(3),R*b1(1),R*b1(2),R*b1(3),0)
quiver3(X2(1),X2(2),X2(3),R*b2(1),R*b2(2),R*b2(3),0)
quiver3(C(1),C(2),C(3),R*b(1),R*b(2),R*b(3),0)
quiver3(C(1),C(2),C(3),R*t1c(1),R*t1c(2),R*t1c(3),0)

daspect([1,1,1])
axis tight



