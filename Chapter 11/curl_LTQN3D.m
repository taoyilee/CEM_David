function curl_LTQN_functions=curl_LTQN3D(lambda,nabla_lambda)
% CURL_LTQN3D returns the twenty LTQN basis functions given lambda and nabla_lamba.
%
% See LTQN3D for further documentation.

global LOCALEDGENODES LOCALFACENODES

curl_LTQN_functions = zeros(20,3);

for iedge = 1:6
    n1 = LOCALEDGENODES(iedge,1);
    n2 = LOCALEDGENODES(iedge,2);
    curl_LTQN_functions(iedge,:)   = 2*cross(nabla_lambda(n1,:),nabla_lambda(n2,:));
    curl_LTQN_functions(iedge+6,:) = 0;
end
for iface = 1:4
    n1 = LOCALFACENODES(iface,1);
    n2 = LOCALFACENODES(iface,2);
    n3 = LOCALFACENODES(iface,3);
    curl_LTQN_functions(iface+12,:) = -3*lambda(n2)*cross(nabla_lambda(n1,:),nabla_lambda(n3,:))...
        -3*lambda(n1)*cross(nabla_lambda(n2,:),nabla_lambda(n3,:));
    curl_LTQN_functions(iface+16,:) = -3*lambda(n3)*cross(nabla_lambda(n2,:),nabla_lambda(n1,:))...
        -3*lambda(n2)*cross(nabla_lambda(n3,:),nabla_lambda(n1,:));
end
