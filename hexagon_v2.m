function vp1=hexagon_v2(a)

spacing = 2*a; 

vp1 = khex(1, 13*spacing, 1); 
vp1 = vp1(2:end, :); %remove the origin

for i = 1:6
    if i < 6
        vp_add = [linspace(vp1(i,1), vp1(i+1,1),14); linspace(vp1(i,2), vp1(i+1, 2), 14)]';
        vp1 = [vp1; vp_add];
    else
        vp_add = [linspace(vp1(i,1), vp1(1,1),14); linspace(vp1(i,2), vp1(1, 2), 14)]';
        vp1 = [vp1; vp_add];
    end
end


vp1 = uniquetol(vp1,'ByRows',true);
