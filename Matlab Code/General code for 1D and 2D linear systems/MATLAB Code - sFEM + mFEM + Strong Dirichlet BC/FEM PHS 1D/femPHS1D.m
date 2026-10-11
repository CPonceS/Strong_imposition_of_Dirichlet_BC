function [model] = femPHS1D(prop,Nelements,TAG_mFEM)

    % TAG_mFEM : interpolation of the strains in the mFEM (Hu-Washizu) model
    %    'C-D' : Continuous displacements - Discontinuous strains (constant strain per element)
    %    'C-C' : Continuous displacements - Continuous strains    (nodal strains, linear)
    
    if nargin < 3 || isempty(TAG_mFEM)
        TAG_mFEM = 'C-D';
    end
    if not(isequal(TAG_mFEM,'C-D') || isequal(TAG_mFEM,'C-C'))
        error('TAG_mFEM must be ''C-D'' or ''C-C''')
    end
    prop.TAG_mFEM = TAG_mFEM;

    % Mesh 1D 
    NE = Nelements;
    NN = NE+1;
    coor  = linspace(prop.dom(1),prop.dom(2),NN); 
    conec = (0:NN-2)'+(1:2);
    Le    = diff(coor)';

    meshFEM.NN    = NN;
    meshFEM.NE    = NE;
    meshFEM.Le    = Le;
    meshFEM.coor  = coor';
    meshFEM.conec = conec;

    % Dimensions and more
    n_PDE   = size(prop.F0,2);    m_PDE   = size(prop.F0,1);
    N_Omega = n_PDE*NN;
    if     isequal(TAG_mFEM,'C-D'); M_Omega = m_PDE*NE;
    elseif isequal(TAG_mFEM,'C-C'); M_Omega = m_PDE*NN;
    end

    Nodes_D = [isequal(prop.Dirichlet,prop.dom(1))*1 isequal(prop.Dirichlet,prop.dom(2))*NN isequal(prop.Dirichlet,[prop.dom(1) prop.dom(2)])*[1 NN] isequal(prop.Dirichlet,[prop.dom(1); prop.dom(2)])*[1 NN]];
    Nodes_D = Nodes_D(Nodes_D>0);
    Nodes_N = [isequal(prop.Neumann,prop.dom(1))*1 isequal(prop.Neumann,prop.dom(2))*NN isequal(prop.Neumann,[prop.dom(1) prop.dom(2)])*[1 NN] isequal(prop.Neumann,[prop.dom(1); prop.dom(2)])*[1 NN]];
    Nodes_N = Nodes_N(Nodes_N>0);
    Nodes_F = setdiff(1:NN,Nodes_D);
    dof_D   = reshape(Nodes_D(:)+(0:n_PDE-1)*NN,[],1);
    dof_N   = reshape(Nodes_N(:)+(0:n_PDE-1)*NN,[],1);
    dof_F   = setdiff(1:n_PDE*NN,dof_D);
    N_free  = n_PDE*(NN-length(Nodes_D));
    N_Dir   = length(dof_D);

    dof_matrix = reshape(1:n_PDE*NN,NN,n_PDE);
    dof_Din = reshape(dof_matrix(Nodes_D,:).*prop.Didx,1,n_PDE*length(Nodes_D)); dof_Din = dof_Din(dof_Din~=0);
    dof_Nin = reshape(dof_matrix(Nodes_N,:).*prop.Nidx,1,n_PDE*length(Nodes_N)); dof_Nin = dof_Nin(dof_Nin~=0);

    [~,idx_sort] = sort([dof_F(:); dof_D(:)]); % sorting for displacement field

    % Input map matrix of Neumann inputs
    [is_free, loc_Nin] = ismember(dof_Nin,dof_F);
    loc_Nin = loc_Nin(is_free); 
    BN = sparse(loc_Nin(:),(1:length(loc_Nin))',1,N_free,length(loc_Nin));

    I_NF  = sparse(eye(N_free));
    I_ND  = sparse(eye(N_Dir));

    NDidx  = length(prop.Didx);  N_bnd  = length(prop.Dirichlet);
    I_full = eye(N_bnd*NDidx);   active_dofs = repelem(logical(prop.Didx),N_bnd);
    Ind_D  = I_full(:,active_dofs);

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
    x1  = coor(conec(:,1))';  
    x2  = coor(conec(:,2))';
    NL  = @(el,x) [(x2(el)-x) (x-x1(el))]/Le(el);
    dNL = [-1 1]./Le;

    Np  = @(el,x) kron(eye(n_PDE),NL(el,x));
    FNp = @(el,x) kron(prop.F0,NL(el,x)) + kron(prop.F1,dNL(el,:));
    if     isequal(TAG_mFEM,'C-D'); Ne = @(el,x) kron(eye(m_PDE),1);
    elseif isequal(TAG_mFEM,'C-C'); Ne = @(el,x) kron(eye(m_PDE),NL(el,x));
    end

    % Allocate memory
    propel = cell(NE,1); 
    % Local matrices            % Global matrices
    M_el  = cell(NE,1);         M  = zeros(N_Omega); 
    Me_el = cell(NE,1);         Me = zeros(M_Omega); 
    Kr_el = cell(NE,1);         Kr = zeros(N_Omega);
    Ke_el = cell(NE,1);         Ke = zeros(M_Omega);
    Fx_el = cell(NE,1);         Fx = zeros(M_Omega,N_Omega);
    
    % Local matrices, vectors and assembly
    NGauss = 2;
    for el = 1:NE
        % Element's properties
        propel{el} = elementProp1D(el,meshFEM,prop);  
        % Local matrices
        M_el{el}  = int_1D( @(x) Np(el,x)'*prop.M(x)*Np(el,x) ,propel{el}.coor(1),propel{el}.coor(end),NGauss);
        Me_el{el} = int_1D( @(x) Ne(el,x)'*Ne(el,x),propel{el}.coor(1),propel{el}.coor(end),NGauss);
        Kr_el{el} = int_1D( @(x) FNp(el,x)'*prop.Ke(x)*FNp(el,x),propel{el}.coor(1),propel{el}.coor(end),NGauss);
        Ke_el{el} = int_1D( @(x) Ne(el,x)'*prop.Ke(x)*Ne(el,x),propel{el}.coor(1),propel{el}.coor(end),NGauss);
        Fx_el{el} = int_1D( @(x) Ne(el,x)'*FNp(el,x),propel{el}.coor(1),propel{el}.coor(end),NGauss);
        % Global matrices
        M  = M  + (propel{el}.Lp'*M_el{el}*propel{el}.Lp);
        Me = Me + (propel{el}.Le'*Me_el{el}*propel{el}.Le);  
        Kr = Kr + (propel{el}.Lp'*Kr_el{el}*propel{el}.Lp);  
        Ke = Ke + (propel{el}.Le'*Ke_el{el}*propel{el}.Le);  
        Fx = Fx + (propel{el}.Le'*Fx_el{el}*propel{el}.Lp);    
    end

    % Partition of matrices 
    Mrr   = sparse(M(dof_F,dof_F));      MrD   = sparse(M(dof_F,dof_D));
    if     isequal(TAG_mFEM,'C-D'); iMe = sparse(diag(1./diag(Me)));   % Me is diagonal
    elseif isequal(TAG_mFEM,'C-C'); iMe = sparse(inv(Me));             % Me is a full (consistent) matrix
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
function propel = elementProp1D(el,meshFEM,prop)

    propel.element  = el; 
    conec_el = meshFEM.conec(el,:);
    coor_el  = meshFEM.coor(conec_el);
    
    % Location matrices (Global to Local)
    n_PDE = size(prop.F0,2); 
    m_PDE = size(prop.F0,1); 
    NN    = length(meshFEM.coor);
    dof_p = n_PDE*NN; 
    Lp = zeros(2*n_PDE,dof_p);    
    Le = []; L0 = zeros(1,NN-1); L0(el) = 1;
    for i=1:n_PDE
        Lp((i-1)*2+1:i*2,(i-1)*NN+((el-1)+1:(el-0)+1)) = eye(2);
    end
    if isequal(prop.TAG_mFEM,'C-D')       % strains discontinuous: one dof per element
        for j=1:m_PDE
            Le = blkdiag(Le,L0);
        end
    elseif isequal(prop.TAG_mFEM,'C-C')   % strains continuous: one dof per node
        Le = zeros(2*m_PDE,m_PDE*NN);
        for j=1:m_PDE
            Le((j-1)*2+1:j*2,(j-1)*NN+(el:el+1)) = eye(2);
        end
    end

    propel.coor  = coor_el;    %  Element node coordinates
    propel.conec = conec_el;   %  Element conectivity
    propel.Lp = Lp;            %  Location matrix for Np(x) interpolated variables
    propel.Le = Le;            %  Location matrix for Ne(x) interpolated variables
end




function I = int_1D(f,a,b,Ngauss)  % integration by Gauss quadrature

    switch Ngauss
        case 1
            w = 2;
            xi= 0;
        case 2
            w =[1 1];
            xi=sqrt(1/3)*[-1 1];
        case 3
            w =[8/9 5/9 5/9];
            xi=[0 sqrt(3/5)*[-1 1]];
        case 4
            w =[((18+sqrt(30))/36)*[1 1] ((18-sqrt(30))/36)*[1 1]];
            xi=[sqrt((3-(2*sqrt(6/5)))/7)*[-1 1] sqrt((3+(2*sqrt(6/5)))/7)*[-1 1]];
        case 5
            w =[128/225 ((322+(13*sqrt(70)))/900)*[1 1] ((322-(13*sqrt(70)))/900)*[1 1]];
            xi=[0 (1/3)*sqrt(5-(2*sqrt(10/7)))*[-1 1] (1/3)*sqrt(5+(2*sqrt(10/7)))*[-1 1]];
        case 10
            w =[0.0666713443*[1 1]  0.1494513492*[1 1]  0.2190863625*[1 1]  0.2692667193*[1 1]  0.2955242247*[1 1]];
            xi=[0.9739065285*[-1 1] 0.8650633667*[-1 1] 0.6794095683*[-1 1] 0.4333953941*[-1 1] 0.1488743390*[-1 1]];
        case 20
            w = [0.152753387130726;0.152753387130726;0.149172986472604;0.149172986472604;0.142096109318382;0.142096109318382;0.131688638449177;0.131688638449177;0.118194531961518;0.118194531961518;0.101930119817240;0.101930119817240;0.0832767415767048;0.0832767415767048;0.0626720483341091;0.0626720483341091;0.0406014298003869;0.0406014298003869;0.0176140071391521;0.0176140071391521];
            xi= [-0.0765265211334973;0.0765265211334973;-0.227785851141645;0.227785851141645;-0.373706088715420;0.373706088715420;-0.510867001950827;0.510867001950827;-0.636053680726515;0.636053680726515;-0.746331906460151;0.746331906460151;-0.839116971822219;0.839116971822219;-0.912234428251326;0.912234428251326;-0.963971927277914;0.963971927277914;-0.993128599185095;0.993128599185095];
    end
    I=0;
    for i=1:Ngauss
        I=I+(w(i)*f(((b-a)*xi(i)+b+a)/2));
    end
    I=((b-a)/2)*I;

end
