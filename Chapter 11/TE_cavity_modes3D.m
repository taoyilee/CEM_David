function [TEeigvalues_sorted,TEmode_sorted] = TE_cavity_modes3D(a,b,d,max_index)
% TE_CAVITY_MODES3D lists the first approximatly TE_z max_index^3 eigenmodes in a rectangular cavity
% See, for instance, Ramo, Whinnery and van Duzer, Chapter 10, "Fields and
% Waves in Communication Electronics",. 3rd edn, Wiley 1994.
index=0;
for mm=0:max_index
    for nn=0:max_index
        for pp=0:max_index
            if any(mm)+any(nn)+any(pp) >= 2% At least two indices must be non-zero
                index=index+1;
                TEeigvalues(index) = sqrt((mm*pi/a)^2+(nn*pi/b)^2+(pp*pi/d)^2);
                TEmode(index,1) = mm;
                TEmode(index,2) = nn;
                TEmode(index,3) = pp;
            end
        end
    end
end
%TEeigvalues;
%TEmode;
[TEeigvalues_sorted,sort_index]=sort(TEeigvalues);
%TEeigvalues_sorted;
TEmode_sorted(:,:)=TEmode(sort_index,:);