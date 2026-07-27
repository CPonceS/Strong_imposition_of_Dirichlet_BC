function meshFEM = gMesh2D(g,Hmax,tag_nodes,tag_edges,tag_elements) % Generate mesh

    mpde = createpde();
    geometryFromEdges(mpde,g);
    order_tag = 'linear'; 
    generateMesh(mpde,'GeometricOrder',order_tag,'Hmin',Hmax,'Hmax',Hmax);
    coor = mpde.Mesh.Nodes'; conec = mpde.Mesh.Elements';
    NE   = size(conec,1);
    
    [~,e,~] = meshToPet(mpde.Mesh);  
    % pdemesh(p,e,t);
    
    pdemesh(mpde, 'NodeLabels', tag_nodes, 'ElementLabels', tag_elements);
    
    numEdges = size(g,2);
    boundaryEdges = cell(1,numEdges);
    for k = 1:numEdges
        idx = find(e(5,:) == k);   % row 5 of 'e' contains the number of edge
        nodes = unique([e(1,idx) e(2,idx)],'stable'); % unique nodes of that edge
        boundaryEdges{k} = nodes;
    end
    
    % Boundary conections
    Bconec = [];
        for i=1:length(boundaryEdges)
            for j=1:length(boundaryEdges{i})-1
                Bconec = [Bconec; [boundaryEdges{i}(j) boundaryEdges{i}(j+1)]];
            end
        end
    
    % Numbering of boundary edges (Global)
    T = zeros(2,size(Bconec,1)); nameE = cell(length(Bconec),1);
    for i=1:length(Bconec)
        T(:,i) = coor(Bconec(i,end),:)'-coor(Bconec(i,1),:)';
        nameE{i,1} = strcat('E',num2str(i));
    end
    if isequal(tag_edges,'on')
        coor_midpoints = T'/2 + coor(Bconec(:,1),:);
        text(coor_midpoints(:,1),coor_midpoints(:,2),nameE,'Color','red')
    end
    hold off
    
    % SHAPE FUNCTIONS
    % coordinates of each node
    x1 = coor(conec(:,1), 1); y1 = coor(conec(:,1), 2);
    x2 = coor(conec(:,2), 1); y2 = coor(conec(:,2), 2);
    x3 = coor(conec(:,3), 1); y3 = coor(conec(:,3), 2);
    
    % Element area
    area = 0.5*abs(x1.*(y2 - y3) + x2.*(y3 - y1) + x3.*(y1 - y2) );
    
    % Coefficients for shape functions
    a = [x2.*y3-x3.*y2,   x3.*y1-x1.*y3,   x1.*y2-x2.*y1];
    b = [y2-y3, y3-y1, y1-y2];
    c = [x3-x2, x1-x3, x2-x1];
    
    % Shape funcions in spatial coordinates (x,y) 
    N    = @(el,x,y) (1/(2*area(el)))*(a(el,:)+b(el,:)*x+c(el,:)*y); 
    dNdx = (1./(2*area)).*b;
    dNdy = (1./(2*area)).*c;
    
    
    % Outputs
    meshFEM.mpde   = mpde;
    meshFEM.g      = g;
    
    meshFEM.ne     = NE; 
    meshFEM.coor   = coor;  
    meshFEM.conec  = conec;
    meshFEM.Bconec = Bconec;
    meshFEM.area   = area;
    
    meshFEM.shapeFcn.N    = N;
    meshFEM.shapeFcn.dNdx = dNdx;
    meshFEM.shapeFcn.dNdy = dNdy;

end