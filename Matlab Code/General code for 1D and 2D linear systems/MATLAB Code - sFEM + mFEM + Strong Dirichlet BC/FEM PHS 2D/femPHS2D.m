function [model] = femPHS2D(prop,meshFEM,TAG_mFEM)

    % TAG_mFEM : interpolation of the strains (Hu-Washizu / mFEM)
    %    'C-C' : Continuous displacements - Continuous strains   (nodal strains, linear)
    %    'C-D' : Continuous displacements - Discontinuous strains (constant strain per element)

    if nargin < 3 || isempty(TAG_mFEM)
        TAG_mFEM = 'C-D';
    end
    if not(isequal(TAG_mFEM,'C-C') || isequal(TAG_mFEM,'C-D'))
        error('TAG_mFEM must be ''C-C'' or ''C-D''')
    end
    prop.TAG_mFEM = TAG_mFEM;    % used by elementPropTriang

    % Dimensions and more
    n_PDE = size(prop.F0,2);     % number of generalized displacements 
    m_PDE = size(prop.F0,1);     % number of generalized strains
    
    NN  = size(meshFEM.coor,1);  % number of nodes
    NE  = size(meshFEM.conec,1); % number of elements
    
    N_Omega = n_PDE*NN;
    if     isequal(TAG_mFEM,'C-D'); M_Omega = m_PDE*NE;
    elseif isequal(TAG_mFEM,'C-C'); M_Omega = m_PDE*NN;
    end
    
    % Boundary conditions
    Nodes_D = unique(meshFEM.Bconec(prop.Dirichlet_Edges,:),'stable');     
    prop.Neumann_Edges = setdiff(1:size(meshFEM.Bconec,1),prop.Dirichlet_Edges')';  % Numbers of global edges where Neumann BC are applied
    nDin = meshFEM.Bconec(prop.Din,:)'; Nodes_Din = unique(nDin(:)','stable');      % Nodes where are applied the Boundary Dirichlet inputs
    nNin = meshFEM.Bconec(prop.Nin,:)'; Nodes_Nin = unique(nNin(:)','stable');      % Nodes where are applied the Boundary Neumann inputs
    
    Nodes_F = setdiff(1:NN,Nodes_D);
    N_free  = n_PDE*(NN-length(Nodes_D));
    
    dof_D = reshape(Nodes_D(:)+(0:n_PDE-1)*NN,[],1);
    dof_F = setdiff(1:N_Omega, dof_D)';
    N_Dir = length(dof_D);
    
    I_NF  = sparse(eye(N_free));
    I_ND  = sparse(eye(N_Dir));
    
    dof_matrix = reshape(1:N_Omega,NN,n_PDE);
    dof_Din = reshape(dof_matrix(Nodes_Din,:).*prop.Didx,1,[]);
    dof_Din = dof_Din(dof_Din ~= 0);
    dof_Nin = reshape(dof_matrix(Nodes_Nin,:).*prop.Nidx,1,[]);
    dof_Nin = dof_Nin(dof_Nin ~= 0);
    
    % Input map matrix of Neumann inputs
    [is_free, loc_Nin] = ismember(dof_Nin,dof_F);
    loc_Nin = loc_Nin(is_free); 
    BN = sparse(loc_Nin(:),(1:length(loc_Nin))',1,N_free,length(loc_Nin));
    
    % Sorting
    Matrix_indD = kron(diag(prop.Didx),eye(length(Nodes_D))); 
    idx_Din = ismember(dof_D,dof_Din);
    Ind_D = sparse(Matrix_indD(:,idx_Din));
    [~,idx_sort] = sort([dof_F(:); dof_D(:)]); % sorting for displacement field
    
    % Useful zero matrices
    Z_NF  = sparse(zeros(N_free));
    Z_FD  = sparse(zeros(N_free,N_Dir));
    Z_DF  = Z_FD';
    Z_DD  = sparse(zeros(N_Dir));
    Z_MM  = sparse(zeros(M_Omega));
    Z_MF  = sparse(zeros(M_Omega,N_free));
    Z_FM  = Z_MF';
    Z_MD  = sparse(zeros(M_Omega,N_Dir));
    Z_DM  = Z_MD';
    
    % Shape functions
    NL     = @(el,x1,x2) meshFEM.shapeFcn.N(el,x1,x2);
    dNLdx1 = meshFEM.shapeFcn.dNdx;
    dNLdx2 = meshFEM.shapeFcn.dNdy;
    
    Np  = @(el,x1,x2) kron(eye(n_PDE),NL(el,x1,x2));
    FNp = @(el,x1,x2) kron(prop.F0,NL(el,x1,x2)) + kron(prop.F1,dNLdx1(el,:)) + kron(prop.F2,dNLdx2(el,:));
    if     isequal(TAG_mFEM,'C-D'); Ne = @(el,x1,x2) kron(eye(m_PDE),1);
    elseif isequal(TAG_mFEM,'C-C'); Ne = @(el,x1,x2) kron(eye(m_PDE),NL(el,x1,x2));
    end
    
    %  Allocate memory
    propel = cell(NE,1); 
    % Local matrices            % Global matrices
    M_el  = cell(NE,1);         M  = zeros(N_Omega); 
    Me_el = cell(NE,1);         Me = zeros(M_Omega); 
    Kr_el = cell(NE,1);         Kr = zeros(N_Omega);
    Ke_el = cell(NE,1);         Ke = zeros(M_Omega);
    Fx_el = cell(NE,1);         Fx = zeros(M_Omega,N_Omega);
    
    % Local matrices, vectors and assembly
    Npoints = 3;
    for el = 1:NE
        % Element's properties
        propel{el} = elementPropTriang(el,meshFEM,prop);  
        % Local matrices
        M_el{el}  = int_triang( @(x1,x2) Np(el,x1,x2)'*prop.M(x1,x2)*Np(el,x1,x2),propel{el}.coor(1:3,1),propel{el}.coor(1:3,2),Npoints);
        Me_el{el} = int_triang( @(x1,x2) Ne(el,x1,x2)'*Ne(el,x1,x2),propel{el}.coor(1:3,1),propel{el}.coor(1:3,2),Npoints);
        Kr_el{el} = int_triang( @(x1,x2) FNp(el,x1,x2)'*prop.Ke(x1,x2)*FNp(el,x1,x2),propel{el}.coor(1:3,1),propel{el}.coor(1:3,2),Npoints);
        Ke_el{el} = int_triang( @(x1,x2) Ne(el,x1,x2)'*prop.Ke(x1,x2)*Ne(el,x1,x2),propel{el}.coor(1:3,1),propel{el}.coor(1:3,2),Npoints);
        Fx_el{el} = int_triang( @(x1,x2) Ne(el,x1,x2)'*FNp(el,x1,x2),propel{el}.coor(1:3,1),propel{el}.coor(1:3,2),Npoints);
        % Global matrices
        M  = M  + (propel{el}.Lp'*M_el{el}*propel{el}.Lp);
        Me = Me + (propel{el}.Le'*Me_el{el}*propel{el}.Le);  
        Kr = Kr + (propel{el}.Lp'*Kr_el{el}*propel{el}.Lp);  
        Ke = Ke + (propel{el}.Le'*Ke_el{el}*propel{el}.Le);  
        Fx = Fx + (propel{el}.Le'*Fx_el{el}*propel{el}.Lp);
    end
    
    % Partition of matrices
    Mrr   = sparse(M(dof_F,dof_F));      MrD = sparse(M(dof_F,dof_D));
    if     isequal(TAG_mFEM,'C-D');      iMe = sparse(diag(1./diag(Me)));   % Me is diagonal
    elseif isequal(TAG_mFEM,'C-C');      iMe = sparse(inv(Me));             % Me is a full (consistent) matrix
    end
    iMrr  = sparse(eye(N_free)/Mrr);      
    Kr_rr = sparse(Kr(dof_F,dof_F));     Kr_rD = sparse(Kr(dof_F,dof_D));
    Kr_Dr = sparse(Kr(dof_D,dof_F));     Kr_DD = sparse(Kr(dof_D,dof_D)); 
    F_mr  = Fx(:,dof_F);                 F_mD  = Fx(:,dof_D);
    
     %% Model sFEM (Hamilton's principle)
    
    model.sFEM.J = sparse([Z_NF -I_NF  Z_FD; 
                           I_NF  Z_NF  Z_FD; 
                           Z_DF  Z_DF  Z_DD]);
    
    model.sFEM.G = sparse([BN                              zeros(N_free,length(dof_D))*Ind_D;
                           zeros(N_free,size(dof_Nin,2))  -iMrr*MrD*Ind_D;
                           zeros(N_Dir ,size(dof_Nin,2))   I_ND*Ind_D]);
    
    model.sFEM.Q = sparse(blkdiag(iMrr,[Kr_rr Kr_rD; Kr_Dr Kr_DD]));
    
    %% Model mFEM (Hu-Washizu principle)
    
    model.mFEM.J = sparse([Z_NF      -F_mr'*iMe' -I_NF    Z_FD; 
                           iMe*F_mr   Z_MM        Z_MF    Z_MD;
                           I_NF       Z_FM        Z_NF    Z_FD  ; 
                           Z_DF       Z_DM        Z_DF    Z_DD]);
    
    model.mFEM.G = sparse([BN                              zeros(N_free,length(dof_D))*Ind_D;
                           zeros(M_Omega,size(dof_Nin,2))  iMe*(F_mD-F_mr*iMrr*MrD)*Ind_D;
                           zeros(N_free,size(dof_Nin,2))  -iMrr*MrD*Ind_D;
                           zeros(N_Dir ,size(dof_Nin,2))   I_ND*Ind_D]);
    
    model.mFEM.Q = sparse(blkdiag(iMrr,Ke,zeros(N_Omega)));
   
    %% Useful outputs
    
    model.TAG_mFEM = TAG_mFEM;
    model.coor     = meshFEM.coor;
    model.NE       = NE;
    model.NN       = NN;
    model.N_free   = N_free;
    model.N_Omega  = N_Omega;
    model.M_Omega  = M_Omega;
    model.idx_sort = idx_sort;

end






%% Auxiliar functions

function propel = elementPropTriang(el,meshFEM,prop) % Properties of the element 'el'

    % Description: Properties of the element 'el'
    %
    % element    %  Number of element
    % coor       %  Element node coordinates
    % conec      %  Element conectivity
    % Lp         %  Location matrix for Np(x,y) interpolated variables
    % Le         %  Location matrix for Ne(x,y) interpolated variables
    
    conec_el = meshFEM.conec(el,:);      
    coor_el  = meshFEM.coor(conec_el,:);
    
    % Location matrices (Global to Local)
    n_PDE = size(prop.F0,2);     m_PDE = size(prop.F0,1);  
    NN = size(meshFEM.coor,1);   NE = size(meshFEM.conec,1); 
    Le = [];  Le0 = zeros(1,NE); Le0(el) = 1;
    % Displacements
    dof_p = n_PDE*NN; 
    Lp = zeros(3*n_PDE,dof_p);
    for i=1:n_PDE
        Lp(3*i-2,(i-1)*NN+conec_el(1)) = 1;
        Lp(3*i-1,(i-1)*NN+conec_el(2)) = 1;
        Lp(3*i-0,(i-1)*NN+conec_el(3)) = 1;
    end  
    if isequal(prop.TAG_mFEM,'C-D')       % strains discontinuous: one dof per element
        for j=1:m_PDE
            Le = blkdiag(Le,Le0);
        end
    elseif isequal(prop.TAG_mFEM,'C-C')   % strains continuous: one dof per node
        Le = zeros(3*m_PDE,m_PDE*NN);
        for j=1:m_PDE
            Le(3*j-2,(j-1)*NN+conec_el(1)) = 1;
            Le(3*j-1,(j-1)*NN+conec_el(2)) = 1;
            Le(3*j-0,(j-1)*NN+conec_el(3)) = 1;
        end
    end
    
    % Definition of the element properties
    propel.element = el;      %  Number of element
    propel.coor  = coor_el;   %  Element node coordinates
    propel.conec = conec_el;  %  Element conectivity
    propel.Lp    = Lp;        %  Location matrix for Np(x,y) interpolated variables
    propel.Le    = Le;        %  Location matrix for Ne(x,y) interpolated variables

end






function I = int_triang(fun,XT,YT,Npoints)  % Integration over triangles

    % Npoints:    Number of points to evaluate the integral
    % fun:        function @(x,y) to be integrated
    % XT:         Column vector with x coordinates of the triangle (counterclockwise)
    % YT:         Column vector with y coordinates of the triangle (counterclockwise)
    
    % I = integral fun(x,y) dxdy
    
    J  = [XT(2)-XT(1) YT(2)-YT(1);
          XT(3)-XT(1) YT(3)-YT(1)];
    
    L3 = @(L1,L2) 1-L1-L2;
    X  = @(L1,L2) [L1 L2 L3(L1,L2)]*XT;
    Y  = @(L1,L2) [L1 L2 L3(L1,L2)]*YT;
    F  = @(L1,L2) fun(X(L1,L2),Y(L1,L2))*det(J);
    
    switch Npoints
        case 1
            l1=1/3;
            l2=1/3;
            w =1;
        case 3
            l1=[1/2 1/2 0];
            l2=[0 1/2 1/2];
            w =[1 1 1]/3;
        case 4
            l1=[1/3 0.6 0.2 0.2];
            l2=[1/3 0.2 0.6 0.2];
            w =[-27 25 25 25]/48;
        case 7
            alpha=[0.797426985353 0.059715871789];
            beta =[0.101286507323 0.470142064105];
            gamma=[0.125939180544 0.132394152788];
            l1=[1/3 alpha(1) beta(1) beta(1) alpha(2) beta(2) beta(2)];
            l2=[1/3 beta(1) alpha(1) beta(1) beta(2) alpha(2) beta(2)];
            w =[0.225 gamma(1) gamma(1) gamma(1) gamma(2) gamma(2) gamma(2)];
    end
    I=0;
    for i=1:length(w)
        I=I+(w(i)*F(l1(i),l2(i)));
    end
    I=0.5*I;
    
    % Reddy - An introduction to the Finite Element Method 3°ED - chap 9, 
    %         table 9.3.2 (pag 559)
end
