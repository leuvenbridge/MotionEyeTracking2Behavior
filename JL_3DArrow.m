function [X,Y,Z] = JL_3DArrow(B,A,r2,r1,L,dots)

% dotted arrow
if nargin<6
   dots=[1 1] ; 
end    
if nargin < 4
   r1=r2*3;L=r2*8; 
end

%
% A = [1 2 0];B=[0 0 1];r1=0.2;r2=0.1;L=0.5;
if L == 0, r1=r2;end


n = dots(1);k=dots(2);if n == 1, k=1;end
l1 = (norm(A-B)-L)/(n-1+k) ; % total length segment+gap ;
l2 = l1*k; % segment length
% [l1 l2 k dots]
[x,y,z]=cylinder([0 r1 r2 r2 0]);
z = repmat([norm(A-B) norm(A-B)-L norm(A-B)-L norm(A-B)-L-l2 norm(A-B)-L-l2]',1,size(z,2));
for i = 1:(n-1)
    [x2,y2,z2]=cylinder([0 r2 r2 0 0]);
    z2 = repmat([0 0 l2 l2 NaN]'+((i-1)*l1),1,size(z,2));
    x=cat(1,x2,x);y=cat(1,y2,y);z=cat(1,z2,z);
end

direction = (A-B)/norm(A-B) ;
rotvec = cross([0 0 1],direction) ;
alpha = acos(dot([0 0 1],direction)) ;
R = vrrotvec2mat([rotvec alpha]) ;
X= x*R(1,1)+y*R(1,2)+z*R(1,3)+B(1);
Y= x*R(2,1)+y*R(2,2)+z*R(2,3)+B(2);
Z= x*R(3,1)+y*R(3,2)+z*R(3,3)+B(3);

% surf(X,Y,Z,'FaceColor',[0 0 0]+0.5)
% xlabel('X');ylabel('Y');zlabel('Z');
% axis equal

