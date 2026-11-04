function [KeB,MeB] = Curved_Beam_3D_Elemental_Stiffness_Matrix( ...
    A,R,Theta,G,Ey,Ez,J,Ix,Iy,Iz,ky,kz,E,rho)
%CURVED_BEAM_3D_ELEMENTAL_STIFFNESS_MATRIX Analytical linear-local arc model.
% REQUIRED INPUTS:
% A: section area; R: centerline radius; Theta: positive arc angle (RADIANS).
% G: shear modulus; Ey,Ez: moduli used for bending about local y,z.
% J: Saint-Venant torsion constant (elastic torsion).
% Ix,Iy,Iz: area moments for rotary inertia; Ix is the POLAR area moment,
%           generally Iy+Iz, and is not generally the torsion constant J.
% OPTIONAL INPUTS (empty also selects default):
% ky=1, kz=1: shear correction factors for y,z displacement.
% E=1: axial Young's modulus; rho=1: mass density.
% Use mutually consistent units. E does NOT multiply Ey,Ez,G again.
%
% Reference curve: r(zeta)=[R*sin(zeta); R*(1-cos(zeta)); 0].
% Local x=tangent, y=inward normal, z=x cross y.
% DOFs at EACH endpoint are expressed in that endpoint's MOVING frame:
% q=[u1 v1 w1 thetaX1 thetaY1 thetaZ1 ...
%    u2 v2 w2 thetaX2 thetaY2 thetaZ2]'.
% Use Transformation_Matrix_Gen_Curved(...,'moving') for global assembly.
%
% N1=1-eta, N2=eta; eta=zeta/Theta, L=R*Theta.
% Strains ordered as [epsilon gammaY gammaZ chiX chiY chiZ]':
% epsilon=u'-v/R; gammaY=v'+u/R-thetaZ; gammaZ=w'+thetaY;
% chiX=thetaX'-thetaY/R; chiY=thetaY'+thetaX/R; chiZ=thetaZ'.
% Prime denotes d/ds. Reference arc is stress free.
%
% KeB: exact integral of B'*D*B for this interpolation and constant section.
% MeB: consistent mass for the SAME independent linear local fields,
%      including axial/transverse inertia and all three rotary inertias.
% No numerical quadrature or Symbolic Math Toolbox is required.
%
% IMPORTANT: This is the BASIC DISPLACEMENT FORMULATION discussed in chat.
% Linear interpolation of moving local components does not reproduce all
% global rigid motions over a finite arc. Exact integration does not fix
% this, nor eliminate locking. Do not treat as a validated production beam.
% A full circle requires closure constraints; these outputs describe an
% open arc with independent endpoint DOFs, including when Theta=2*pi.
% Does not include geometric stiffness, pressure, thermal loading or warping.

L = R*Theta;

%% Section constitutive matrix
D = diag([E*A, ky*G*A, kz*G*A, G*J, Ey*Iy, Ez*Iz]);

%% Separate B(eta) = C + eta*H
% Build full B at eta=0 and eta=1. The difference is the eta coefficient.
Ns = [-1,1]/L;
Bstart = local_B([1,0],Ns,R);
Bend   = local_B([0,1],Ns,R);
C = Bstart;
H = Bend-Bstart;

%% One-point (midpoint) integration -- removes membrane/shear locking
Bmid = C + 0.5*H;
KeB  = L*(Bmid.'*D*Bmid);
KeB  = 0.5*(KeB+KeB.');

%% Consistent mass: integral of interpolation'*inertia*interpolation
% Inertia is diagonal for centroidal principal section axes.
% Integral N1^2=Integral N2^2=1/3, Integral N1*N2=1/6.
Inertia = rho*diag([A,A,A,Ix,Iy,Iz]);
MeB = (L/6)*[2*Inertia, Inertia; Inertia, 2*Inertia];
end

function B = local_B(N,Ns,R)
B = zeros(6,12);
for i = 1:2
    Bi = [Ns(i), -N(i)/R, 0,     0,      0,       0; ...
          N(i)/R, Ns(i), 0,     0,      0,      -N(i); ...
          0,      0,     Ns(i), 0,      N(i),    0; ...
          0,      0,     0,     Ns(i), -N(i)/R, 0; ...
          0,      0,     0,     N(i)/R, Ns(i),   0; ...
          0,      0,     0,     0,      0,       Ns(i)];
    idx = (i-1)*6+(1:6);
    B(:,idx) = Bi;
end
end
