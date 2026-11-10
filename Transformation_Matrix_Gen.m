function [T,a] = Transformation_Matrix_Gen(X1,X2,vref)

if nargin < 3 || isempty(vref)
  vref = [0 0 1];
end


d = X2(:)' - X1(:)';
L = norm(d);

if L == 0
  error('Transformation_Matrix_Gen: the two nodes coincide.');
end

a = L/2;

ex = d/L;                                  % local x: along the beam

% local z: reference direction with its axial part removed
ez = vref(:)' - dot(vref(:)',ex)*ex;
if norm(ez) < 1e-8*norm(vref)              % beam parallel to vref: fall back
    vref = [0 1 0];
    ez = vref - dot(vref,ex)*ex;
end
ez = ez/norm(ez);

ey = cross(ez,ex);                         % local y completes right-hand set

R = [ex; ey; ez];                          % rows = local axes in global coords
T = blkdiag(R,R,R,R);                      % u, theta at node 1 and node 2
end
