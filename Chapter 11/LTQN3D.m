function LTQN_functions=LTQN3D(lambda,nabla_lambda)
% LTQN3D returns the twenty LTQN basis functions given lambda and nabla_lambda.
%
% The first six entries are the Whitney functions,
% the next six the LTLN functions, and the last eight
% the LTQN "face" functions.
% Note that the associations
% edge 1: 12
% edge 2: 13
% edge 3: 14
% edge 4: 23
% edge 5: 24
% edge 6: 34
%
% face 1: 123
% face 2: 124
% face 3: 134
% face 4: 234
% are hard-wired into this function via the global variables.

global LOCALEDGENODES LOCALFACENODES

LTQN_functions = zeros(20,3);
for iedge = 1:6
    n1 = LOCALEDGENODES(iedge,1);
    n2 = LOCALEDGENODES(iedge,2);
    LTQN_functions(iedge,:)    = lambda(n1)*nabla_lambda(n2,:) - lambda(n2)*nabla_lambda(n1,:);
    LTQN_functions(iedge+6,:)  = lambda(n1)*nabla_lambda(n2,:) + lambda(n2)*nabla_lambda(n1,:);
end
for iface = 1:4
    n1 = LOCALFACENODES(iface,1);
    n2 = LOCALFACENODES(iface,2);
    n3 = LOCALFACENODES(iface,3);
    LTQN_functions(iface+12,:) = lambda(n2)*lambda(n3)*nabla_lambda(n1,:)+ lambda(n1)*lambda(n3)*nabla_lambda(n2,:)...
        -2*lambda(n1)*lambda(n2)*nabla_lambda(n3,:);
    LTQN_functions(iface+16,:) = lambda(n3)*lambda(n1)*nabla_lambda(n2,:)+ lambda(n2)*lambda(n1)*nabla_lambda(n3,:)...
        -2*lambda(n2)*lambda(n3)*nabla_lambda(n1,:);
end
