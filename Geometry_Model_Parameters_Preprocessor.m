function[ndof,nel,neq,nnp,LM]=Geometry_Model_Parameters_Preprocessor(Nodes,IEN,Type)

% fprintf('   Called Function: Geometry_Model_Parameters_Preprocessor\n')   
% fprintf('          Location:%s\n',pwd)
% fprintf('          Purpose: Set basic FEM model properties including LM matrix,\n\n')


%% Initialize System  
    nsd     = 3;                    % number of space dimensions 
    ndof 	= 6;                    % number of degrees-of-freedom per node
     
    nnp 	= Nodes(end,1);     	% number of nodal points
    S=size(IEN);
    nel 	= S(2);                 % number of elements
    neq     = ndof*nnp;             % number of equations    
    
%% generate LM array 
    LM=zeros(max(Type(2,:))*ndof,nel);
    for e = 1:nel
        nen=Type(2,e);
        for j = 1:nen
            for m = 1:ndof
                ind = (j-1)*ndof + m;
                 LM(ind,e) = ndof*IEN(j,e) - ndof + m;
            end
        end
    end   

%% Generate Predefined Model Matrix
    % [IEN_Colm,IEN_Colm2]=Predefined_Matrix_Assembly(nel,Type,IEN);
    
    % fprintf('Number of Nodes=%d\nNumer of Elements=%d\nNumer of Equations=%d\n',nnp,nel,neq)

end

function[IEN_Colm,IEN_Colm2]=Predefined_Matrix_Assembly(nel,Type,IEN)

    %% IEN_Colm is column representation of Element nodes, used for 1D variable such as Temperature
        IEN_Colm=zeros(sum(Type(2,:)),1);
        count=1;
        for e=1:nel
            IEN_Colm(count:count+Type(2,e)-1,1)=IEN(1:Type(2,e),e);
            count=count+Type(2,e);
        end
    
    %% IEN_Colm2 is column representation of Element nodes, used for multi-variable such as Displacement [u1;v1]
        IEN_Colm2=zeros(length(IEN_Colm)*3,1);       % IEN_Stress
        count=1;
        for i=1:length(IEN_Colm)
            IEN_Colm2(count)=3*IEN_Colm(i)-2;
            IEN_Colm2(count+1)=3*IEN_Colm(i)-1;
            IEN_Colm2(count+2)=3*IEN_Colm(i);
            count=count+3;
        end
        
%     %% Assembly_F maps IEN_Colm2 from local to Global and Assemble them, like adding Force for node 3 of element 1 and node 3 for element 3 together
%         Assembly_F=zeros(neq,length(IEN_Colm2));    
%         for i=1:length(IEN_Colm2)
%             Assembly_F(IEN_Colm2(i),i)=1;
%         end 
%         
%     %% Assembly_S maps IEN_Colm from local to Global and Assemble them, like adding Stress for node 3 of element 1 and node 3 for element 3 together
%         Assembly_S=zeros(nnp,length(IEN_Colm));    
%         for i=1:length(IEN_Colm)
%             Assembly_S(IEN_Colm(i),i)=1;
%         end 
        
%     %%
%     Elem_sh=zeros(nnp,1);
%     for i=1:nnp
%         Elem_sh(i,1)=sum(Assembly_S(i,:));
%     end
%         
end
        
        