function ypr = JL_R2yapirod(R)

if size(R,3)>1
   ypr = zeros(size(R,1),3);
   for i = 1:size(R,1)
       ypr(i,:) = JL_R2yapirod(squeeze(R(i,:,:)));
   end
else
    % rotation matrix from euler angles, in the sequence
    % first about x-axis, then y-axis, finally z-axis
    % R = Rz(y)Ry(p)Rx(r)
    % R = [cos(p)*cos(y) cos(y)*sin(p)*sin(r)-cos(r)*sin(y) cos(r)*cos(y)*sin(p)+sin(r)*sin(y); ...
    %      cos(p)*sin(y) sin(p)*sin(r)*sin(y)+cos(r)*cos(y) cos(r)*sin(p)*sin(y)-cos(y)*sin(r); ...
    %        -sin(p)              cos(p)*sin(r)                     cos(p)*cos(r)];
    
    y = atan2d(R(2,1),R(1,1));
    p = -asind(R(3,1));
    r = atan2d(R(3,2),R(3,3));
%     r = yapiro(y,p,0);
%     x = r\R;
%     r = atan2d(x(3,2),x(2,2));

    ypr = [y p r];
end