function meshFEM = gMesh2D_aux(g,Hmax,tag_nodes,tag_edges,tag_elements)

    numEdges = size(g,2);

    % Mesh
    P = []; edgeID = [];  segs = cell(1,numEdges); 
    nOff = 0;
    for k = 1:numEdges
        x1 = g(2,k); x2 = g(3,k); y1 = g(4,k); y2 = g(5,k);
        if g(1,k) == 2                       % straight line
            L = hypot(x2-x1,y2-y1);
            n = max(1,ceil(L/Hmax));
            t = linspace(0,1,n+1)';
            pts = [x1+t*(x2-x1), y1+t*(y2-y1)];
        elseif g(1,k) == 1 || g(1,k) == 4    % circle arc o elipse
            xc = g(8,k); yc = g(9,k);
            if g(1,k) == 1
                ea = g(10,k); eb = g(10,k); phi = 0;
            else
                ea = g(10,k); eb = g(11,k); phi = g(12,k);
            end
            Rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];
            % p = c + Rot*[ea*cos(t); eb*sin(t)]
            d0 = Rot'*[x1-xc; y1-yc];  t0 = atan2(d0(2)/eb, d0(1)/ea);
            d1 = Rot'*[x2-xc; y2-yc];  t1 = atan2(d1(2)/eb, d1(1)/ea);
            if t1 <= t0, t1 = t1 + 2*pi; end
            ts  = linspace(t0,t1,200);
            xy  = Rot*[ea*cos(ts); eb*sin(ts)];
            Larc = sum(hypot(diff(xy(1,:)),diff(xy(2,:))));
            n  = max(2,ceil(Larc/Hmax));
            th = linspace(t0,t1,n+1);
            xy = Rot*[ea*cos(th); eb*sin(th)];
            pts = [xc+xy(1,:)', yc+xy(2,:)'];
            pts(1,:)   = [x1 y1];            
            pts(end,:) = [x2 y2];          
        else
            error('Type of edge %d not supported',g(1,k));
        end
        P = [P; pts];
        segs{k} = (nOff+1:nOff+size(pts,1))';
        nOff = nOff + size(pts,1);
    end

    [Pu,~,ic] = uniquetol(P,1e-10,'ByRows',true);
    boundaryEdges = cell(1,numEdges);
    C = [];
    for k = 1:numEdges
        ids = ic(segs{k})';
        boundaryEdges{k} = ids;
        C = [C; ids(1:end-1)' ids(2:end)'];
    end

    % Triangulation
    DT0 = delaunayTriangulation(Pu,C);
    inside0 = isInterior(DT0);

    xmin = min(Pu(:,1)); xmax = max(Pu(:,1));
    ymin = min(Pu(:,2)); ymax = max(Pu(:,2));
    dy = Hmax*sqrt(3)/2;
    Q = [];
    for j = 0:ceil((ymax-ymin)/dy)
        off = mod(j,2)*Hmax/2;
        xs = (xmin+off:Hmax:xmax)';
        Q = [Q; xs, (ymin+j*dy)*ones(size(xs))];
    end
    tri = pointLocation(DT0,Q);
    ok  = ~isnan(tri);
    ok(ok) = inside0(tri(ok));                       
    Q = Q(ok,:);
        DTp = delaunayTriangulation(Pu);             
    vid = nearestNeighbor(DTp,Q);                    
    far = hypot(Q(:,1)-DTp.Points(vid,1), Q(:,2)-DTp.Points(vid,2)) > 0.7*Hmax;
    Q = Q(far,:);

    % Final mesh
    DT = delaunayTriangulation([Pu;Q],C);
    conec = DT.ConnectivityList(isInterior(DT),:);
    coor  = DT.Points;
    NE    = size(conec,1);

    x = @(i) coor(conec(:,i),1); y = @(i) coor(conec(:,i),2);
    A2 = (x(2)-x(1)).*(y(3)-y(1)) - (x(3)-x(1)).*(y(2)-y(1));
    neg = A2 < 0;
    conec(neg,[2 3]) = conec(neg,[3 2]);

    % Plot
    figure('color','w')
    triplot(conec,coor(:,1),coor(:,2),'color',[0.278 0.6 0.812]); axis equal; hold on
    if isequal(tag_nodes,'on')
        text(coor(:,1),coor(:,2),strcat('n',string(1:size(coor,1))),'Color','k')
    end
    if isequal(tag_elements,'on')
        cen = (coor(conec(:,1),:)+coor(conec(:,2),:)+coor(conec(:,3),:))/3;
        text(cen(:,1),cen(:,2),strcat('e',string(1:NE)),'color',[0.2 0.55 0.2],'HorizontalAlignment','center')
    end

    % Conectivity at the boundary
    Bconec = [];
    for i = 1:length(boundaryEdges)
        Bconec = [Bconec; boundaryEdges{i}(1:end-1)' boundaryEdges{i}(2:end)'];
    end
    % Boundary plot
    xb = [coor(Bconec(:,1),1) coor(Bconec(:,2),1) nan(size(Bconec,1),1)]';
    yb = [coor(Bconec(:,1),2) coor(Bconec(:,2),2) nan(size(Bconec,1),1)]';
    plot(xb(:), yb(:), 'color', [0.9 0.5 0.2], 'LineWidth', 0.5);

    T = zeros(2,size(Bconec,1)); nameE = cell(size(Bconec,1),1);
    for i = 1:size(Bconec,1)
        T(:,i) = coor(Bconec(i,end),:)' - coor(Bconec(i,1),:)';
        nameE{i,1} = strcat('E',num2str(i));
    end
    if isequal(tag_edges,'on')
        coor_midpoints = T'/2 + coor(Bconec(:,1),:);
        text(coor_midpoints(:,1),coor_midpoints(:,2),nameE,'color','red')
    end
    hold off

    % Boundary plot handle
    xb = [coor(Bconec(:,1),1) coor(Bconec(:,2),1) nan(size(Bconec,1),1)]';
    yb = [coor(Bconec(:,1),2) coor(Bconec(:,2),2) nan(size(Bconec,1),1)]';
    Bplot = @() plotBoundary(xb(:),yb(:));

    % Shape functions
    x1 = coor(conec(:,1),1); y1 = coor(conec(:,1),2);
    x2 = coor(conec(:,2),1); y2 = coor(conec(:,2),2);
    x3 = coor(conec(:,3),1); y3 = coor(conec(:,3),2);
    area = 0.5*abs(x1.*(y2-y3) + x2.*(y3-y1) + x3.*(y1-y2));
    a = [x2.*y3-x3.*y2, x3.*y1-x1.*y3, x1.*y2-x2.*y1];
    b = [y2-y3, y3-y1, y1-y2];
    c = [x3-x2, x1-x3, x2-x1];
    N    = @(el,x,y) (1/(2*area(el)))*(a(el,:)+b(el,:)*x+c(el,:)*y);
    dNdx = (1./(2*area)).*b;
    dNdy = (1./(2*area)).*c;

    % Outputs
    meshFEM.g      = g;
    meshFEM.ne     = NE;
    meshFEM.coor   = coor;
    meshFEM.conec  = conec;
    meshFEM.Bconec = Bconec;
    meshFEM.Bplot  = Bplot;
    meshFEM.area   = area;
    meshFEM.shapeFcn.N    = N;
    meshFEM.shapeFcn.dNdx = dNdx;
    meshFEM.shapeFcn.dNdy = dNdy;
end


% Auxiliar function
function plotBoundary(xb,yb)

    wasHeld = ishold;
    hold on
    plot3(xb, yb, zeros(numel(xb),1), 'color', [1 1 1]*0.6, 'LineWidth', 0.5);
    if ~wasHeld, hold off; end
end