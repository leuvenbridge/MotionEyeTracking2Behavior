function Az = JL_R2TAz(R)    

v = cross(R(:,3),[0 0 1]');
if R(9)> 0
    Rf = (vrrotvec2mat([v; real(asin(norm(v)))]))*R(:,1);
else
    Rf = (vrrotvec2mat([v; pi-real(asin(norm(v)))]))*R(:,1); % rotate the HD vector according to tilt
end
Az = atan2d(Rf(2),Rf(1));