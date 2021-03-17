function [h] = Monkey_3D_V2(x,m,R)
if (nargin <1), x = [1 1 1] ; end
if (nargin <2), m = [0 0 0] ; end
if (nargin <3), R = eye(3,3) ; end
if length(x) ==1 ; x = [x x x]; end

[xs,ys,zs] = sphere(35) ;

% axis equal; grid on
d = 0.2 ;
hold_state = ishold(gca) ;
hold on
% clear h
h(1) = surf(xs,ys,zs,36*ones(size(xs)),'EdgeColor','none','AlphaData',0,'FaceColor',[0.5 0.5 0.5]) ;

[xs,ys,zs] = sphere(12) ;
% one eye
h(2) = surf(xs*d +0.8,ys*d -0.34,zs*d +0.35,60*ones(size(xs)),'EdgeColor','none','FaceColor',[1 1 1]*0.9) ;
h(3) = surf(xs*d/2+d*0.8+0.8,ys*d/2-0.34,zs*d/2 +0.35,0*ones(size(xs)),'EdgeColor','none','FaceColor',[0.1 0.1 0.1]) ;

% other eye
h(4) = surf(xs*d +0.8,ys*d +0.34,zs*d +0.35,64*ones(size(xs)),'EdgeColor','none','FaceColor',[1 1 1]*0.9) ;
h(5) = surf(xs*d/2+d*0.8+0.8,ys*d/2+0.34,zs*d/2 +0.35,0*ones(size(xs)),'EdgeColor','none','FaceColor',[0.1 0.1 0.1]) ;

[xs,ys,zs] = sphere(35) ;

xm = xs*0.4; ym = ys*0.7 ; zm = zs*0.5  ;
a = -15 ;
% mouth
tmp = xm*cosd(a)-zm*sind(a) ; zm = zm*cosd(a)+xm*sind(a) ; xm = tmp ;
xm = xm+0.75 ; zm = zm-0.4 ;
h(6) = surf(xm,ym,zm,55*ones(size(xs)),'EdgeColor','none','FaceColor',[1 1 1]*0.8) ;

xm = xs*0.05 +0.95 ; ym = ys*0.57 ; zm = zs*0.4 -0.42  ;
% h(9) = surf(xm,ym,zm,55*ones(size(xs)),'EdgeColor','none','FaceColor',[0.1 0.1 0.1]) ;

% hair
xm = xs*0.99 -0.07 ; ym = ys*0.95 ; zm = zs*0.95 + 0.13  ;
h(7) = surf(xm,ym,zm,55*ones(size(xs)),'EdgeColor','none','FaceColor',[1 1 1]*0.2) ;

xm = xs*0.98 - 0.1; ym = ys*0.98 ; zm = zs*0.9 + 0.05  ;
h(8) = surf(xm,ym,zm,55*ones(size(xs)),'EdgeColor','none','FaceColor',[1 1 1]*0.2) ;
% set(gca,'CameraPosition',[17 -3 4])

% [z_,y_,x_] = JL_3DRotArrow(0,27,0.8,0.03,0.03,0,1);
% h(end+1)=surf(x_+1.1,y_,-z_+0.1,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
% h(end+1)=surf(x_+1.1,-y_,-z_+0.1,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
% dx=1.1;dy=mean(y_(end,:)); dz=-mean(z_(end,:))+0.1;
% [x_,y_,z_]=sphere(20);x_=x_*0.03+dx;y_=y_*0.03+dy;z_=z_*0.03+dz;
% h(end+1)=surf(x_,y_,z_,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
% h(end+1)=surf(x_,-y_,z_,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
%   
% JL_RM_Objects(h(9:end),yapirod(0,14,0),[0.16 0 0.27]) ;

[z_,y_,x_] = JL_3DRotArrow(0,34,0.8,0.03,0.03,0,1);z_=z_*0.8;
h(end+1)=surf(x_+1.1,y_,-z_+0.1,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
h(end+1)=surf(x_+1.1,-y_,-z_+0.1,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
dx=1.1;dy=mean(y_(end,:)); dz=-mean(z_(end,:))+0.1;
[x_,y_,z_]=sphere(20);x_=x_*0.03+dx;y_=y_*0.03+dy;z_=z_*0.03+dz;
h(end+1)=surf(x_,y_,z_,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
h(end+1)=surf(x_,-y_,z_,'EdgeColor','none','FaceColor',[0 0 0]+0.2) ;
  
JL_RM_Objects(h(9:end),yapirod(0,-45,0),[-0.06 0 -1.05]) ;



% for i = 1:length(h)
%    set(h(i),'XData', get(h(i),'XData')*x(1)+m(1)) ;
%    set(h(i),'YData', get(h(i),'YData')*x(2)+m(2)) ;
%    set(h(i),'ZData', get(h(i),'ZData')*x(3)+m(3)) ;
%    set(h(i),'CDataMapping','direct') ;
%
% end
if not(hold_state), hold off;end

for i = 1:length(h)
    
    X = get(h(i),'XData')*x(1) ;
    Y = get(h(i),'YData')*x(2) ;
    Z = get(h(i),'ZData')*x(3) ;
    set(h(i),'XData', X*R(1,1)+Y*R(1,2)+Z*R(1,3)+m(1)) ;
    set(h(i),'YData', X*R(2,1)+Y*R(2,2)+Z*R(2,3)+m(2)) ;
    set(h(i),'ZData', X*R(3,1)+Y*R(3,2)+Z*R(3,3)+m(3)) ;
    set(h(i),'CDataMapping','direct') ;
    
end
