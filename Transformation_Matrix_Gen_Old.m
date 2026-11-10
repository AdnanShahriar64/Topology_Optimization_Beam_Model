function[T,a]=Transformation_Matrix_Gen_Old(Nodes,IEN,Frame_Or,e)

% if e==1
%     fprintf('       Called Function (several times): Transformation_Matrix_Gen \n')
%     fprintf('              Location:%s\n',pwd)
%     fprintf('              Purpose: Obtain Transformation Matrix for each Elements\n\n')
% end

x2=Nodes(IEN(2,e),2);
x1=Nodes(IEN(1,e),2);

y2=Nodes(IEN(2,e),3);
y1=Nodes(IEN(1,e),3);

z2=Nodes(IEN(2,e),4);
z1=Nodes(IEN(1,e),4);

L=((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^0.5;
a=L/2;

%% Orientation

    Vy1=[x1 y1 z1];
    Vy2=[x2 y2 z2];

    Vx=Vy2-Vy1;
    Vxu=Vx/norm(Vx);

    if Frame_Or(1)==0       % Radial direction is y
        Vz=cross(Vx,Vy1);
        Vzu=Vz/norm(Vz);
        % fprintf('Considering y axis along radial axis')
    elseif Frame_Or(1)==1   % system in x y plane
        Vz=[0 0 1];
        Vzu=Vz/norm(Vz);
    end

    % Vyu=cross(Vxu,Vzu);

    lx=Vxu(1);
    mx=Vxu(2);
    nx=Vxu(3);

    lz=Vzu(1);
    mz=Vzu(2);
    nz=Vzu(3);

    ly=mz*nx-nz*mx;
    my=nz*lx-lz*nx;
    ny=lz*mx-mz*lx;

    T3=[lx mx nx ; ly my ny; lz mz nz];

    T=zeros(12,12);
    for i=1:4
        T(i*3-2:i*3,i*3-2:i*3)=T3;
    end

end
