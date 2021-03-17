function behavior = fmmExtractBehavior(deuteron, spike2)

% NOTE: the timestamps here are not synced with spike timestamps yet
ts_markers  = spike2.headpos.ts;
ts_delta = 0.02;
ts       = 0:ts_delta:ts_markers(end);  % resample at 50 Hz for all behavioral variables

expdate  = str2double(deuteron.expt.date);

% unit: mm
% Attention!!! center and front markers were swapped due
% to setting in Vicon during calibration, keep it this way then
nMarker = 4;
nDim = 3;
markers_raw  = cat(3, spike2.headpos.center.xyz', ...  % front marker
                      spike2.headpos.front.xyz',...    % back marker
                      spike2.headpos.left.xyz',...     % left marker
                      spike2.headpos.right.xyz');      % right marker                   

markers_raw = shiftdim(markers_raw, 2);
nSample = size(markers_raw, 3);
ok_missing = markers_raw(1,1,:).^2==0 | markers_raw(2,1,:).^2==0 | markers_raw(3,1,:).^2==0 | markers_raw(4,1,:).^2==0;    % where data is missing
markers_raw(:,:,ok_missing) = nan;
%
miss_st_all = [];
miss_end_all = [];
for iMarker = 1:nMarker
    nok = isnan(squeeze(markers_raw(iMarker,1,:)));
    miss_st = find([0;diff(nok)]==1)-2;  % add 2 extract samples before and after the gap
    miss_end = find([diff(nok);0]==-1)+2;
    if ~isempty(miss_st) && ~isempty(miss_end)
        % make sure the start and end of missing trunk are paired
        if miss_end(1)-2 < miss_st(1)+2
            miss_end(1) = [];            
        end
        if miss_st(end)+2 > miss_end(end)-2
            miss_st(end) = [];
        end
    end
    miss_st_all = cat(1, miss_st_all, miss_st);
    miss_end_all = cat(1, miss_end_all, miss_end);
    for i = 1:length(miss_st)
        markers_raw(:,:,miss_st(i):miss_end(i)) = nan;
    end
end
markers_raw = markers_raw(:,:,1:nSample);
%

markers = fillmissing(markers_raw, 'spline', 3, 'EndValues', 'none');

%
% detect gaps that are larger than 0.5 s and replace with nan 
% because these gaps can't be filled correctly
for i = 1:length(miss_st_all)
    if miss_end_all(i)-miss_st_all(i) > 100  % 1 s
        markers(:,:,miss_st_all(i):miss_end_all(i)) = nan;
    end
end
%

% interpolate to the new ts
markers_intp = nan(nMarker, nDim, length(ts));
for i = 1:nMarker
    for j = 1:nDim
        markers_intp(i, j, :) = interp1(ts_markers, squeeze(markers(i,j,:)), ts);
    end
end

%%
h_arena = 2120;   % height of the arena, mm
r_arena = 1650;   % radius of the arena, mm

% translation matrix from marker center to head center
% head center is in the same egocentric horizontal plane as the eye center
if strcmp(deuteron.expt.monkey, 'Kraut')
    trans_markerhead = [0, 0, -126.0];        % head center relative to the marker center, in head coordinate, measured from MRI & CT images
    trans_markereye  = [50.2, 17.7, -126.0];  % eye center relative to the marker center, in head coordinate, measured from MRI & CT images
elseif strcmp(deuteron.expt.monkey, 'Lysander')
    trans_markerhead = [0, 0, -106.0];
    trans_markereye  = [nan, nan, nan];
elseif strcmp(deuteron.expt.monkey, 'Bruno')
    trans_markerhead = [0, 0, -105.1];
    trans_markereye  = [49.4, 20.5, -105.1];
end

% extract the rotation matrix from arena to head reference frames
% translate marker center to head center and eye center, respectively
% i.e. head reference frames with origin at head center and eye center
markers_headcenter = nan(size(markers_intp));
markers_eyecenter  = nan(size(markers_intp));
% rotation matrix from arena to head coodinate
R_arenatohead = nan(3,3,size(markers_intp,3));
% gravity vector in head reference frame
gravity_head  = nan(3,size(markers_intp,3));
% tilted azimuth angle using JL code
azimuth_tilt  = nan(1,size(markers_intp,3));
% yaw, pitch, roll velocity using JL code
yawpitchroll_vel = nan(3,size(markers_intp,3)) ;

% head ori. & arena intersection
head_arena_inter = nan(3, size(markers_headcenter,3));
% coordinate in the unfolded 2D arena map
head_arena_2dpos = nan(2, size(markers_headcenter,3));

%%
for i = 1:size(markers_intp,3)
    % head center
    m = squeeze(markers_intp(:,:,i));
    % get the z axis vector in head reference frame
    vv = cross(m(1,:)-m(3,:),m(1,:)-m(4,:));
    if vv(3) < 0     % v1(3) should be always negative, i.e. vector pointing upword
        % sometimes Vicon mistakes left and right markers
        % because of the symmetry
        markers_intp([3 4],:,i) = markers_intp([4 3],:,i);
    end
    m = squeeze(markers_intp(:,:,i));

    if sum(isnan(m(:)))==0
    %%
        % ex, ey, and ez define the head reference frame
        % ex: back to front; ey: right to left; ez: upward
        ey = m(3,:)-m(4,:); ey = ey./norm(ey);
        v1 = m(1,:)-m(3,:); v2 = m(1,:)-m(4,:);
        ez = cross(v1, v2); ez = ez./norm(ez);
        ex = cross(ey,ez);  ex = ex./norm(ex);
        
        % this is the rotation matrix at this moment
        r   = [ex; ey; ez];
        R_arenatohead(:,:,i) = r;
        % gravity vector in head
        gh  = r * [0;0;-1];
        gravity_head(:,i) = gh./norm(gh);
        azimuth_tilt(i)   = JL_R2TAz(r');

        % ****** NOTE FROM JL
        % I would have computed directly 
        % markers_headcenter(i,:) = mean(m) + (r.'*trans_markerhead')';
        % i.e. the 3D coordinates of head center
        % that seems more intuitive
        % but in the end it give the same result as [x1 y1 z1] below

        % translate the markers to head center and eye center:
        % convert the translation vector (from markers' center to head and eye center,
        % in head reference frame, see above) to arena reference frame and add the
        % resulting vector to markers coordinate matrix
        markers_headcenter(:,:,i) = m + (r.'*trans_markerhead')';
        markers_eyecenter(:,:,i)  = m + (r.'*trans_markereye')';
        % rotation velocity
        if i > 1
            % given rotation matrix at timestamp i (Ri) and at timestamp i-1 (R(i-1))
            % the rotation matrix from Ri to R(i-1) is: (Ri)'*R(i-1)  (Ra->b = Ra->o * Ro->b = (Ro->a)'*Ro->b)
            % inv(R) = R' since the rows/columns are orthogonal
            if isnan(R_arenatohead(1,1,i)) || isnan(R_arenatohead(1,1,i-1))
                yawpitchroll_vel(:,i) = [nan; nan; nan];
            else
                % inv(R_arenatohead(:,:,i)') = R_arenatohead(:,:,i); inv is inaccurate
%                 yawpitchroll_vel(:,i) = (-JL_R2yapirod(inv(R_arenatohead(:,:,i)')*R_arenatohead(:,:,i-1)')/ts_delta)';
                yawpitchroll_vel(:,i) = (-JL_R2yapirod(R_arenatohead(:,:,i)'\R_arenatohead(:,:,i-1)')/ts_delta)';  % it's the same results
            end
        end

        % where the head points at on the inner surface of the cylinder arena
        % i.e. where the head orientation vector intersects with the arena
        % use this to approximate spatial view locations for data without eye tracking

        markers_headcenter_mean = squeeze(nanmean(markers_headcenter(:,:,i),1));
        x1 = markers_headcenter_mean(1);
        y1 = markers_headcenter_mean(2);
        z1 = markers_headcenter_mean(3);
        % head orientation vector in arena reference frame
        hv = ex;

        % ****** Note from JL: all this is almost good. 
        % Btw you can use https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection if you have a more complex plane equation
        % also Matlab will solve polynomials for you ;) (function roots) 
        % finally you can also write the line with a parameter k 
        % (x, y, z) = (xl, yl, zl) + k*hv
        % and write the polynom x.^2+y.^2 = r_arena.^2 as a second order
        % polynom with variable k. Then solve the polynom and select solution
        % with k>0

        % head center to head-arena intersection (x,y,z) vector
        % is parallel to the head orientation vector hv
        % the intersection meets either x.^2+y.^2 = r_arena.^2 or z = 0 or z = h_arena
        % (x-x1)/(y-y1) = hv(1)/hv(2), substitute y in x.^2+y.^2 = r_arena.^2
        % and get it in this form: a*x.^2 + b*x + c = 0
        a = 1+(hv(2)/hv(1)).^2;
        b = 2*hv(2)/hv(1)*(y1-hv(2)/hv(1)*x1);
        c = (y1-hv(2)/hv(1)*x1).^2-r_arena.^2;
        rs = b.^2 - 4*a*c;

        % if there is real solution
        if rs >= 0
            % analytical solutions
            xval = [(-b+sqrt(rs))/2/a, (-b-sqrt(rs))/2/a];
            % only one solution corresponds to where the head points
            % the other solution corresponds to intersection behind the animal
            if (x1-xval(1))*(-hv(1)) > 0
                xx = xval(1);
            else
                xx = xval(2);
            end
            % solve z using: (x-x1)/(z-z1) = hv(1)/hv(3)
            zz = hv(3)/hv(1)*(xx-x1)+z1;
            % when the intersection is below the floor or above the ceiling
            if zz > h_arena
                zz = h_arena;
            elseif zz < 0
                zz = 0;
            end
            % redo x and y using: (x-x1)/(z-z1) = hv(1)/hv(3) and
            % (y-y1)/(z-z1) = hv(2)/hv(3)
            xx = hv(1)/hv(3)*(zz-z1)+x1;
            yy = hv(2)/hv(3)*(zz-z1)+y1;
            head_arena_inter(:,i) = [xx; yy; zz];

            % view position in the unfolded 2D arena map
            % half height at 0 degree is the origin:
            %                           %%%%%
            %                          %     %       # Ceiling
            %                          %     %
            %                           %   %
            %            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %            %                                  %
            %            %                O                 %    # Wall
            %            %                                  %
            %            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           %   %
            %                          %     %
            %                          %     %       # floor
            %                           %%%%%

            if zz == h_arena
                head_arena_2dpos(1,i) = xx; 
                head_arena_2dpos(2,i) = 1/2 * h_arena + r_arena - yy;
            elseif zz == 0
                head_arena_2dpos(1,i) = xx; 
                head_arena_2dpos(2,i) = -1/2 * h_arena - r_arena + yy;
            else
                if xx >= 0
                    head_arena_2dpos(1,i) = acos(yy/sqrt(xx.^2+yy.^2)) * r_arena;
                else
                    head_arena_2dpos(1,i) = -acos(yy/sqrt(xx.^2+yy.^2)) * r_arena;
                end
                head_arena_2dpos(2,i) = zz - 1/2 * h_arena;
            end
        end

        if (0)
            % This section will draw a beautiful plot for you to check
            % Look for the necessary functions in a folder
            % (JL_Standard_Matlab_Scripts) on the lab's server
            subplot(121) ;
            cla
            hold on
            hh = Monkey_3D_V2(100) ;
            jll = 3000 ;
            [jlx, jly,jlz]=JL_3DArrow([0 0 0],[jll 0 0],25,0,0,[10 0.5]);
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 0 0],'EdgeColor','none') ;
            JL_RM_Objects(hh,r',[x1 y1 z1]) ;

            [jlx, jly,jlz]=cylinder([r_arena, r_arena],50) ;jlz=jlz*h_arena ;
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 0 0]+0.75,'EdgeColor','none') ;

            jll = 1000 ;
            [jlx, jly,jlz]=JL_3DArrow([0 0 0],[jll 0 0],25,80,120);
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 0 1],'EdgeColor','none') ;
            hh(end+1) = surf(jly,jlx,jlz,'FaceColor',[0 0.5 0],'EdgeColor','none') ;
            hh(end+1) = surf(jlz,jly,jlx,'FaceColor',[1 0 0],'EdgeColor','none') ;

            [jlx, jly,jlz]=sphere(30);
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 1 1],'EdgeColor','none') ;
            JL_RM_Objects(hh(end),100,[xx yy zz]) ;


            light('Position',[100 -100 300]*1000)
            set(hh,'AmbientStrength',0.7,'DiffuseStrength',1,'SpecularStrength',0.3,'FaceLighting','gouraud') ;
            axis equal

            xlabel('X')
            ylabel('Y')
            zlabel('Z')

            subplot(224) ;
            cla
            hold on
            plot(head_arena_2dpos(1,:),head_arena_2dpos(2,:),'.k')
            plot(head_arena_2dpos(1,i),head_arena_2dpos(2,i),'or','MarkerFaceColor','r')
            axis equal
        end
    end

end

%% allocentric position and elevation in arena
position    = squeeze(nanmean(markers_headcenter, 1));
position_xy = position(1:2,:);
elevation_z = position(3,:);

% translation speed in 3D
pos_diff    = diff(position,1,2);
dist_diff   = sqrt(pos_diff(1,:).^2+pos_diff(2,:).^2+pos_diff(3,:).^2);
trans_speed = [nan dist_diff./ts_delta];
trans_speed_sm = medfilt1(trans_speed,5);
% movement direction in earth horizontal
mv_dir_earth = [nan atan2d(pos_diff(2,:), pos_diff(1,:))];

% head direction vector
head_dir_vector = squeeze(markers_headcenter(1,:,:) - markers_headcenter(2,:,:));
% azimuth head direction in earth horizontal
azimuth_earth = atan2d(head_dir_vector(2,:), head_dir_vector(1,:));

% arena center to head center vector
u = position_xy;
% head direction vector in horizontal
v = head_dir_vector(1:2,:);

% egocentric boundary angle, [-180, 180] degrees
% angle between u = [x1,y1] and v = [x2,y2] equals atan2d(x1*y2-y1*x2,x1*x2+y1*y2)
% measures in ccw from u to v [0 180], or cw [-180 0]
a_eb = -atan2d(u(1,:).*v(2,:)-u(2,:).*v(1,:), u(1,:).*v(1,:)+u(2,:).*v(2,:));
% egocentric boundary distance
% radius minus arena center to head center distance
d_eb = r_arena - sqrt(u(1,:).^2+u(2,:).^2);
% polar to cartesian coordinate
[x_eb,y_eb] = pol2cart(a_eb/180*pi, d_eb);

%% output structure
behavior.ts                     = ts;
behavior.position_xy            = position_xy;
behavior.elevation_z            = elevation_z;
behavior.trans_speed            = trans_speed;
behavior.trans_speed_sm         = trans_speed_sm;
behavior.headdir_azimuth_tilt   = azimuth_tilt;
behavior.movedir_azimuth_earth  = mv_dir_earth;
behavior.headdir_azimuth_earth  = azimuth_earth;
behavior.tilt                   = gravity_head;
behavior.yawpitchroll_vel       = yawpitchroll_vel;
behavior.head_arena_inter       = head_arena_inter;
behavior.head_arena_2dpos       = head_arena_2dpos;
behavior.ego_boundary_polar     = [a_eb; d_eb];
behavior.ego_boundary_xy        = [x_eb; y_eb];

%% eye data
% no eye data for Lysander
% eye data for Kraut and Bruno
if isfield(spike2, 'eye') && std(spike2.eye.lhor) > 0.1  % check if there is eye tracking data
    % eye direction in horizontal plane
    eyelhor = spike2.eye.lhor;
    eyelver = spike2.eye.lver;
    % date when switched to mirror-based eye tracker, horizontal fliped, vertical the same
    if expdate >= 20190906      
        eyelhor = -eyelhor;
    end
    % some preprocessing of eye data
    % pay attention to + / -
    if isfield(spike2.eye, 'deltat')
        dt_eye = spike2.eye.deltat;
    else
%         dt_eye = 0.001;   % 1k Hz sampling rate set in spike2 script, actual sampling rate could fluctuate
        dt_eye = ts_markers(end)/length(eyelhor);
    end
    % filter eye data by thresholding eye movemnet speed to remove potential artifact
    deyelhor = diff(eyelhor);
    nok = abs(deyelhor) > 20;
    % take one extra point before and after the nan gaps
    nok(2:end) = nok(2:end) | nok(1:end-1);
    nok(1:end-1) = nok(1:end-1) | nok(2:end);
    nok = [false;nok];
    eyelhor(nok) = nan;
    eyelver(nok) = nan;
    % remove data outside +- 40 degree in both directions
    % potential artifact
    eyelhor(abs(eyelhor)>40) = nan;
    eyelver(abs(eyelver)>40) = nan;

    ts_eye = 0:dt_eye:length(eyelhor)*dt_eye-dt_eye;
    % sometimes spike2 outputs one fewer datapoint in either hor or ver
    if length(eyelhor) < length(ts_eye)
        eyelhor = [eyelhor;eyelhor(end)];
    end
    if length(eyelver) < length(ts_eye)
        eyelver = [eyelver;eyelver(end)];
    end
    % resample at 50 Hz
    ex = interp1(ts_eye, eyelhor, ts);
    ey = interp1(ts_eye, eyelver, ts);
    
    % eye orientation vector in head
    e_v_head = [cosd(ex); sind(ex); tand(ey)];

    eyedir_azimuth_earth = nan(1, size(markers_eyecenter,3));
    eye_arena_inter = nan(3, size(markers_eyecenter,3));
    eye_arena_2dpos = nan(2, size(markers_eyecenter,3));
    markers_eyecenter_mean = squeeze(nanmean(markers_eyecenter,1));

    for i = 1:size(markers_eyecenter,3)

        R_this = squeeze(R_arenatohead(:,:,i));
        % eye orientation vector in arena reference frame
        ev = R_this.' * e_v_head(:,i);
        eyedir_azimuth_earth(i) = atan2d(ev(2,:), ev(1,:));
        
        % solve eye-arena intersection the same way as
        % head-arena intersection shown above
        % eye center to eye-arena intersection vector
        % is parallel to the eye orientation vector
        % the intersection meets either x.^2+y.^2=r_arena.^2 or z = 0 or h_arena
        % (x-x1)/(y-y1) = ev(1)/ev(2)
        x1 = markers_eyecenter_mean(1,i);
        y1 = markers_eyecenter_mean(2,i);
        z1 = markers_eyecenter_mean(3,i);
        a = 1+(ev(2)/ev(1)).^2;
        b = 2*ev(2)/ev(1)*(y1-ev(2)/ev(1)*x1);
        c = (y1-ev(2)/ev(1)*x1).^2-r_arena.^2;
        rs = b.^2 - 4*a*c;

        % if there is real solution
        if rs >= 0
            xval = [(-b+sqrt(rs))/2/a, (-b-sqrt(rs))/2/a];
            if (x1-xval(1))*(-ev(1)) > 0
                xx = xval(1);
            else
                xx = xval(2);
            end
            zz = ev(3)/ev(1)*(xx-x1)+z1;
            if zz > h_arena
                zz = h_arena;
            elseif zz < 0
                zz = 0;
            end
            xx = ev(1)/ev(3)*(zz-z1)+x1;
            yy = ev(2)/ev(3)*(zz-z1)+y1;
            eye_arena_inter(:,i) = [xx; yy; zz];
            
            % view position in the unfolded 2D arena map
            % half height at 0 degree is the origin
            if zz == h_arena
                eye_arena_2dpos(1,i) = xx; 
                eye_arena_2dpos(2,i) = 1/2 * h_arena + r_arena - yy;
            elseif zz == 0
                eye_arena_2dpos(1,i) = xx; 
                eye_arena_2dpos(2,i) = -1/2 * h_arena - r_arena + yy;
            else
                if xx >= 0
                    eye_arena_2dpos(1,i) = acos(yy/sqrt(xx.^2+yy.^2)) * r_arena;
                else
                    eye_arena_2dpos(1,i) = -acos(yy/sqrt(xx.^2+yy.^2)) * r_arena;
                end
                eye_arena_2dpos(2,i) = zz - 1/2 * h_arena;
            end
        end
        
        if (0) % Change to if(1) to draw 
            % This section will draw a beautiful plot for you to check
            % Look for the necessary functions in a folder
            % (JL_Standard_Matlab_Scripts) on the lab's server
            % Note that the monkey's head size is eggagerated. The center of
            % the eye corresponds to your data but the center of the head will
            % be a bit off....
            
            subplot(121) ;
            cla
            hold on
            hh = Monkey_3D_V2(100) ;
            JL_RM_Objects(hh,1,[-nanmean([hh(4).XData(:);hh(5).XData(:)]) -nanmean([hh(4).YData(:);hh(5).YData(:)]) -nanmean([hh(4).ZData(:);hh(5).ZData(:)])]) ;
            jll = 3000 ;
            [jlx, jly,jlz]=JL_3DArrow([0 0 0],e_v_head(:,i)'*jll,25,0,0,[10 0.5]);
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 0 0],'EdgeColor','none') ;
            JL_RM_Objects(hh,R_this.',[x1 y1 z1]) ;
            
            [jlx, jly,jlz]=cylinder([r_arena, r_arena],50) ;jlz=jlz*h_arena ;
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 0 0]+0.75,'EdgeColor','none') ;
            
            jll = 1000 ;
            [jlx, jly,jlz]=JL_3DArrow([0 0 0],[jll 0 0],25,80,120);
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 0 1],'EdgeColor','none') ;
            hh(end+1) = surf(jly,jlx,jlz,'FaceColor',[0 0.5 0],'EdgeColor','none') ;
            hh(end+1) = surf(jlz,jly,jlx,'FaceColor',[1 0 0],'EdgeColor','none') ;
            
            [jlx, jly,jlz]=sphere(30);
            hh(end+1) = surf(jlx,jly,jlz,'FaceColor',[0 1 1],'EdgeColor','none') ;
            JL_RM_Objects(hh(end),100,[xx yy zz]) ;
            
            
            light('Position',[100 -100 300]*1000)
            set(hh,'AmbientStrength',0.7,'DiffuseStrength',1,'SpecularStrength',0.3,'FaceLighting','gouraud') ;
            axis equal
            
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            
            subplot(222) ;
            cla
            hold on
            plot(ex,ey,'.k')
            plot(ex(i),ey(i),'or','MarkerFaceColor','r')
            axis equal
            xlabel('eye hor (°)');ylabel('eye vert (°)')
            subplot(224) ;
            cla
            hold on
            plot(eye_arena_2dpos(1,:),eye_arena_2dpos(2,:),'.k')
            plot(eye_arena_2dpos(1,i),eye_arena_2dpos(2,i),'or','MarkerFaceColor','r')
            axis equal
            
            pause(0.01)
        end
        
    end
    
    behavior.eyepos_xy              = [ex; ey];    % eye in head
    behavior.eyedir_azimuth_earth   = eyedir_azimuth_earth;
    behavior.eye_arena_inter        = eye_arena_inter;
    behavior.eye_arena_2dpos        = eye_arena_2dpos;
end

end