function[Nodes_L,IEN_L,normal_plane_vector] = Discretized_Local_Geomtry_for_Curved_Beam_Element(e, lead, Nodes, NB, nArc)

  ndofe = 12;

  # Extract internal nodes
  X1 = Nodes(NB(2,e),2:4);
  X2 = Nodes(NB(3,e),2:4);

  # Element chord
  c  = X2 - X1;
  Lc = norm(c);
  ch = c/Lc;

  # Correcting lead
  if lead < 0
      lead = -1
  elseif lead > 0
      lead =1
  end


  # --- end tangents (zero if neighbor missing) ---
  t1 = [0 0 0];  t2 = [0 0 0];
  if NB(1,e) > 0
    X0 = Nodes(NB(1,e),2:4);
    t1 = (X1 - X0)/norm(X1 - X0);
  end
  if NB(4,e) > 0
    X3 = Nodes(NB(4,e),2:4);
    t2 = (X3 - X2)/norm(X3 - X2);
  end

  # --- weights from lead (lead=0 -> symmetric least squares) ---
  w1 = (1 - lead)/2;
  w2 = (1 + lead)/2;

  if norm(t1) == 0
    w1 = 0; w2 = 1;
  end
  if norm(t2) == 0, w2 = 0
    w1 = 1;
  end


  # --- optimal plane: n = perpendicular part of d ---
  d     = w1*t1 - w2*t2;
  dperp = d - dot(d,ch)*ch;
  B     = norm(dperp);

  if B < 1e-12
      error('Tangents along the chord: plane undefined.');
  end

  n = dperp/B;                    # in-plane normal to chord (old t1c)
  b = cross(n, ch);               # bending-plane normal (old b)

  # --- optimal half-sweep angle ---
  A     = dot(w1*t1 + w2*t2, ch);
  alpha = atan2(B, A);            # 0..pi

  a = Lc/2;
  M = (X1 + X2)/2;
  if sin(alpha) < 1e-8
      error('Optimal arc is straight (alpha ~ 0): use a straight element.');
  end
  R = a/sin(alpha);
  h = -a*cot(alpha);
  C = M + h*n;                    # circle center on the bisector

  # --- report misfit ---
  T1 = cos(alpha)*ch + sin(alpha)*n;
  T2 = cos(alpha)*ch - sin(alpha)*n;
  if w1 > 0, fprintf('entry mismatch = %f deg\n', acosd(min(1,dot(T1,t1)))); end
  if w2 > 0, fprintf('exit  mismatch = %f deg\n', acosd(min(1,dot(T2,t2)))); end

  # --- nodes: psi from -alpha (X1) to +alpha (X2) ---
  Nodes_L = zeros(nArc+1,4);
  Nodes_L(1,:) = [1 X1];
  Nodes_L(2,:) = [2 X2];
  for i = 1:nArc-1
      psi = -alpha + 2*alpha*i/nArc;
      Nodes_L(i+2,1)   = i+2;
      Nodes_L(i+2,2:4) = C + R*(sin(psi)*ch + cos(psi)*n);
  end



  IEN_L = zeros(2,nArc);
  IEN_L(1,1) = 1;
  for i = 1:nArc-1
      IEN_L(2,i)   = i+2;
      IEN_L(1,i+1) = i+2;
  end
  IEN_L(2,nArc) = 2;



  normal_plane_vector = b;



  %figure
  %ord = [1 3:nArc+1 2];
  %plot3(Nodes_L(ord,2),Nodes_L(ord,3),Nodes_L(ord,4))
  %hold on
  %scatter3(Nodes_L(:,2),Nodes_L(:,3),Nodes_L(:,4))
  %try
  %  plot3([X0(1) X1(1)],[X0(2) X1(2)],[X0(3) X1(3)],'k')
  %  quiver3(X1(1),X1(2),X1(3),R*b1(1),R*b1(2),R*b1(3),0)
  %end
  %try
  %  plot3([X2(1) X3(1)],[X2(2) X3(2)],[X2(3) X3(3)],'k')
  %  quiver3(X2(1),X2(2),X2(3),R*b2(1),R*b2(2),R*b2(3),0)
  %end


  %quiver3(C(1),C(2),C(3),R*b(1),R*b(2),R*b(3),0)
  %quiver3(C(1),C(2),C(3),R*n(1),R*n(2),R*n(3),0)

  %daspect([1,1,1])
  %axis tight


end
