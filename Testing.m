clear
clc

Nodes = [1  0     0     0;
         2  0     0.1  -1;
         3  0.5  -0.1  -1.1
         4  0.5  -0.1   0.1];





Nodes = [1  0     0     0;
         2  0     0.1  -1;
         3  0.5   0.1  -1.5;
         4  1     0.1  -1.5;
         5  1.25  0.1  -1.25;
         6  1.25  0.1   0];


[IEN, Type, beta] = Generate_IEN_Type(Nodes);

sadasd

IEN = [Nodes(1:end-1,1)'; Nodes(2:end,1)'];

Type = zeros(4,size(IEN,2));
Type(1,:)=1;
Type(1,2)=2;
Type(1,4)=2;

nel = size(IEN,2);


Geometry.Circular.Outer_Radius = 0.009525;
Geometry.Circular.Thickness    = 0.0008890;


Mat_Prop.MoE   = 120e9;
Mat_Prop.PoRat = 0.34;
Mat_Prop.rho   = 8940;


NB      = zeros(4,nel);
for e = 1:nel

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


% [Ke,Me] = Condensed_Beam_Elemental_Matrix_Generator(2,Nodes,NB,Type,Geometry,Mat_Prop);


lead = 0;
nArc = 16;


figure
hold on
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),'filled')

ord = [1 3:nArc+1 2];

for e = 1:size(IEN,2)
    if Type(1,e) == 1                      % straight: just connect the two nodes
        P = Nodes(IEN(:,e),2:4);
        plot3(P(:,1),P(:,2),P(:,3),'k','LineWidth',1.5)
    else                                   % curved: discretized arc
        [Nodes_L,IEN_L,normal_plane_vector] = Discretized_Local_Geomtry_for_Curved_Beam_Element(e, lead, Nodes, NB, nArc);



        plot3(Nodes_L(ord,2),Nodes_L(ord,3),Nodes_L(ord,4),'b','LineWidth',1.5)
        scatter3(Nodes_L(:,2),Nodes_L(:,3),Nodes_L(:,4),10,'b')
    end
end

daspect([1 1 1])
axis tight
view(3)

