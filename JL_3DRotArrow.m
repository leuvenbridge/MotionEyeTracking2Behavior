function [X,Y,Z] = JL_3DRotArrow(a1,a2,r,r1,r2,L,direction)


%%
a1 = a1*sign(direction);a2=a2*sign(direction);
% da = mod(a2-a1,360) ;% r = 2 ;
da=a2-a1;
% r2 = 0.2;
% r1 = 0.1;
% L = 0.5/2 ;
% a1 = 190 ; a2 = -30;


[x,z,y]=cylinder([0 r1*ones(1,da)]);x=x+r;y=y*0;
   
y(1,:)=x(1,:)*sind(-da) ;
x(1,:)=x(1,:)*cosd(-da) ;
for i = 2:da+1
   y(i,:)=x(i,:)*sind(i-da-2) ;
   x(i,:)=x(i,:)*cosd(i-da-2) ;    
end
% 
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         if x(i,j)>0&y(i,j)>-L&y(i,j)<0;
%             x(i,j) = sqrt(x(i,j)^2+y(i,j)^2-L^2) ;
%             y(i,j)=-L;
%         end
%     end
% end
I = find(x>0&y>-L&y<0) ;
x(I)=sqrt(x(I).^2+y(I).^2-L^2) ;
y(I)=-L;

[x2,z2,y2]=cylinder([r2 0]) ;x2=x2+r;
y2(1,:)=-L ;y2(2,:)=L;
x = [x;x2];y=[y;y2];z=[z;z2];

y = y*sign(direction);
a1 = a1*sign(direction);a2=a2*sign(direction);

X = x*cosd(a2)-y*sind(a2);
Y = y*cosd(a2)+x*sind(a2);
Z=z;

% surf(X,Y,Z,'FaceColor',[0 0 0]+0.5)
% xlabel('X');ylabel('Y');zlabel('Z');
% axis equal
%%
% A = [1 2 0];B=[0 0 1];r1=0.2;r2=0.1;L=0.5;
% [x,y,z]=cylinder([0 r1 r2 r2 0]);z = repmat([norm(A-B) norm(A-B)-L norm(A-B)-L 0 0]',1,size(z,2));
% 
% direction = (A-B)/norm(A-B) ;
% rotvec = cross([0 0 1],direction) ;
% alpha = acos(dot([0 0 1],direction)) ;
% rotvec = rotvec/norm(rotvec) ;
% R = vrrotvec2mat([rotvec alpha]) ;
% X= x*R(1,1)+y*R(1,2)+z*R(1,3)+B(1);
% Y= x*R(2,1)+y*R(2,2)+z*R(2,3)+B(2);
% Z= x*R(3,1)+y*R(3,2)+z*R(3,3)+B(3);

% surf(X,Y,Z,'FaceColor',[0 0 0]+0.5)
% xlabel('X');ylabel('Y');zlabel('Z');
% axis equal

