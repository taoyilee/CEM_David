function facemake3D
% FACEMAKE Make faces on mesh.
% Element 1 first:
global ELEMENTS ELEMENT_FACES NUM_ELEMENTS FACES NUM_FACES LOCALFACENODES
FACES(1,:) = [ELEMENTS(1,1),ELEMENTS(1,2),ELEMENTS(1,3)];
FACES(2,:) = [ELEMENTS(1,1),ELEMENTS(1,2),ELEMENTS(1,4)];
FACES(3,:) = [ELEMENTS(1,1),ELEMENTS(1,3),ELEMENTS(1,4)];
FACES(4,:) = [ELEMENTS(1,2),ELEMENTS(1,3),ELEMENTS(1,4)];
ELEMENT_FACES(1,1) = 1;
ELEMENT_FACES(1,2) = 2;
ELEMENT_FACES(1,3) = 3;
ELEMENT_FACES(1,4) = 4;
% Now other elements
face_counter = 4;
for ielem = 2:NUM_ELEMENTS
   for jface = 1:4 % a tet has four faces
     TEMPFACES(1,:) = [ELEMENTS(ielem,LOCALFACENODES(jface,1)),ELEMENTS(ielem,LOCALFACENODES(jface,2)),ELEMENTS(ielem,LOCALFACENODES(jface,3))];
     new_face = 1; % Default: true. Re-set for each face.
     % Now test if this face has already been assigned
     for kface = 1:face_counter % check all previously assigned faces
       if (TEMPFACES(1,:) == FACES(kface,:)) 
         new_face = 0;
         ELEMENT_FACES(ielem,jface) = kface;
         break;
       end
     end
     if new_face
       face_counter = face_counter+1;
       FACES(face_counter,:) = TEMPFACES(1,:);
       ELEMENT_FACES(ielem,jface) = face_counter;
     end
   end
end
NUM_FACES = face_counter;
