function [K] = Assembly_Truss(Ke,Type,ndof,LM,nel,neq)


%     fprintf('       Called Function: Assembly_Truss\n')
%     fprintf('              Location:%s\n',pwd)
%     fprintf('              Purpose: Obtain Global Stiffness or Mass matrix for Truss or Beam section\n\n')




K 	= zeros(neq);

%% Assembly of K and M
%fprintf('\nAssemblying...')

for e = 1:nel
    for loop1 = 1:Type(2,e)*ndof
        i = LM(loop1,e);
        for loop2 =  1:Type(2,e)*ndof
            j = LM(loop2,e);
            K(i,j) = K(i,j) + Ke(loop1,loop2,e);
        end
    end
    %fprintf('Element=%d   Type=%d\n',e,Type(e));
end
%fprintf('Done\n')


end
