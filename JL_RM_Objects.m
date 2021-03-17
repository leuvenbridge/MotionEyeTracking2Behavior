function [h] = JL_RM_Objects(h,R,M)
if nargin < 3, M = [0 0 0]; end
if numel(R)==1, R = eye(3,3)*R;end
if numel(R)==3, R = yapirod(R);end
for i = 1:length(h)
    x = get(h(i),'XData');
    y = get(h(i),'YData');
    z = get(h(i),'ZData');
    X= x*R(1,1)+y*R(1,2)+z*R(1,3)+M(1);
    Y= x*R(2,1)+y*R(2,2)+z*R(2,3)+M(2);
    Z= x*R(3,1)+y*R(3,2)+z*R(3,3)+M(3);
    
    set(h(i),'XData',X);
    set(h(i),'YData',Y);
    set(h(i),'ZData',Z);
end