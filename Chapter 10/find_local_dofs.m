function find_local_dofs(dof_RWG))
%FIND_LOCAL_DOFS finds the two faces per edge-based dof and
% saves them in datastructure DOFLOCALNUM 

global  ELEMENT_EDGES NUM_ELEMENTS NUM_EDGES DOFLOCALNUM 

for iedge = 1:NUM_EDGES
    edge_found = 0; % flag
    terminate = 0;
    for jelem = 1:NUM_ELEMENTS
        for kedge = 1:3
            if ELEMENT_EDGES(jelem,kedge) == iedge
                if ~edge_found
                    DOFLOCALNUM(dof_RWG(iedge),1) = kedge;
                    edge_found = 1;
                else
                    DOFLOCALNUM(dof_RWG(iedge),2) = kedge;
                    terminate = 1; % Terminate element loop; max of two faces per edge
                    break
                end
            end
        end
        if terminate
            break % exit element loop
        end
    end
end

end

