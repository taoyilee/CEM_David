function plot_field(dofs,dof_e1,XX,YY,rows,cols,plotnum,eigvalue)
% PLOT_FIELD plots the vector field at the grid points defined
% vectors XX and YY.

global ELEMENTS NODE_COORD EDGES ELEMENT_EDGES NUM_ELEMENTS LOCALEDGENODES NUM_DOFS
for ii = 1:length(XX)
    for jj = 1:length(YY)
      point_found = 0; 
      for i_elem = 1:NUM_ELEMENTS,
        trinodes = ELEMENTS(i_elem,:); 
        x1 = NODE_COORD(trinodes(1),1);
        y1 = NODE_COORD(trinodes(1),2);
        x2 = NODE_COORD(trinodes(2),1);
        y2 = NODE_COORD(trinodes(2),2);
        x3 = NODE_COORD(trinodes(3),1);
        y3 = NODE_COORD(trinodes(3),2);
        % Find bounding rectangle
        xlower = min([x1 x2 x3]);
        xupper = max([x1 x2 x3]);
        ylower = min([y1 y2 y3]);
        yupper = max([y1 y2 y3]);
        % Get coordinates for this point 
        xc = XX(ii);
        yc = YY(jj);
        % Check if point lies inside bounding rectangle
        if (xlower <= xc & xc <= xupper) & (ylower <= yc & yc <= yupper)
          % if is does, then compute simplex coordinates and check their validity.   
          lambda=simplex2D(i_elem,xc,yc);
          tolerance = 1e-4;
          if -tolerance <= lambda(1) & lambda(1) <= 1+tolerance & ...
             -tolerance <= lambda(2) & lambda(2) <= 1+tolerance & ...
             -tolerance <= lambda(3) & lambda(3) <= 1+tolerance
            % This point lies in or on this element. 
            edges = ELEMENT_EDGES(i_elem,1:3);
            % Get global dof's associated with this triangle.
            % if 0, then prescribed
            dof_e1_tri(1) = dof_e1(edges(1)); 
            dof_e1_tri(2) = dof_e1(edges(2)); 
            dof_e1_tri(3) = dof_e1(edges(3)); 
            %keyboard
            % Now get the global dof associated with the local edge.
            for k_edge = 1:3,
              if dof_e1_tri(k_edge)
                % keyboard
                dofs_tri(k_edge) = dofs(dof_e1_tri(k_edge));
              else
                 dofs_tri(k_edge) = 0;  % prescribed
              end,
            end
            ctln_funcs = whitney(i_elem,xc,yc);
            vec_field(ii,jj,:) =  dofs_tri(1)*ctln_funcs(1,:) + ...
                                  dofs_tri(2)*ctln_funcs(2,:) + ...
                                  dofs_tri(3)*ctln_funcs(3,:); 
            
            %keyboard;
            point_found = 1;
            break
          end
		end
      end
      if ~point_found
          keyboard % problem here
      end  
    end
end
[Xgrid,Ygrid] = meshgrid(XX,YY);
% Note: the transposes below are needed due to the convention used by the function QUIVER. 
subplot(rows,cols,plotnum),quiver(Xgrid',Ygrid',vec_field(:,:,1),vec_field(:,:,2));
title(['k_c =',num2str(eigvalue)])