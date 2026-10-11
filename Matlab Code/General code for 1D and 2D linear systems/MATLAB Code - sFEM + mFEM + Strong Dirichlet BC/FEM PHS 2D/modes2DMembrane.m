function modes2DMembrane(model,meshFEM,method,N_modes,View,Boundary_ref,FAC)
%
% Plots the first N_modes eigenmodes of a 2D Membrane
%
% The surface is colored by the displacement magnitude
%
%   model        : output of femPHS2D
%   meshFEM      : output of gMesh2D  (coor, conec, Bplot)
%   method       : 'sFEM' or 'mFEM'
%   N_modes      : number of modes to be plotted
%   View         : [Az El]                  (optional, default [35 30])
%   Boundary_ref : 'on' or 'off'            (optional, default 'on')
%   FAC          : amplifying factor        (optional, default 1)
%
% State vector (its dimension depends on the model):
%   sFEM:  x = [ p (N_free) | r_free (N_free) | r_D (N_Dir) ]
%   mFEM:  x = [ p (N_free) | eps (M_Omega) | r_free (N_free) | r_D (N_Dir) ]
%
% The displacements are the last N_Omega entries, ordered as [r_free; r_D].
% model.idx_sort puts them in nodal order: r = [psi_1 (NN); psi_2 (NN); w (NN)].

    if nargin < 5 || isempty(View);          View = [35 30];      end
    if nargin < 6 || isempty(Boundary_ref);  Boundary_ref = 'on'; end
    if nargin < 7 || isempty(FAC);           FAC = 1;             end

    A = full(model.(method).J*model.(method).Q);

    NN      = model.NN;
    N_free  = model.N_free;
    N_Omega = model.N_Omega;
    M_Omega = model.M_Omega;
    n_PDE   = N_Omega/NN;
    coor    = meshFEM.coor;
    conec   = meshFEM.conec;

    if     isequal(method,'mFEM'); factor = 1;
    elseif isequal(method,'sFEM'); factor = 0;
    else;  error('method must be ''sFEM'' or ''mFEM''')
    end

    idx_r = N_free + factor*M_Omega + 1:size(A,1);   % displacements [r_free; r_D]

    % Eigenvalues and eigenvectors
    [V,lam] = eig(A);
    f       = imag(diag(lam))/(2*pi);        % [Hz]
    f_min   = 0.001;                         % discards zero eigenvalues
    idx     = find(f > f_min);               % one eigenvalue of each conjugate pair
    [f,order] = sort(f(idx));  idx = idx(order);
    N_modes = min(N_modes,numel(f));

    w_max = (max(coor(:,2))-min(coor(:,2)))/4;     % plotted amplitude
    cmap  = custom_colormap(256,0.165);

    % Single figure: square or rectangle with the fewest rows that fits N_modes
    nc = ceil(sqrt(N_modes));                % columns
    nr = ceil(N_modes/nc);                   % rows
    figure('Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.8]);

    for mode = 1:N_modes

        r = V(idx_r,idx(mode));
        w = r(model.idx_sort);                         % transverse displacement

        [~,imax] = max(abs(w));
        w = real(w*exp(-1i*angle(w(imax))));           % align the phase -> real mode shape
        coorZ = FAC*(w/max(abs(w)))*w_max;

        ax = subplot(nr,nc,mode);                      % row-wise: left to right, top to bottom
        hold(ax,'on'); box(ax,'on')
        trisurf(conec,coor(:,1),coor(:,2),coorZ,coorZ,'Parent',ax);   % Eigenmode
        if isequal(Boundary_ref,'on') && isfield(meshFEM,'Bplot')
            meshFEM.Bplot();                           % Boundary reference (z = 0)
        end
        hold(ax,'off')

        TITLE = strcat({'Mode '},num2str(mode),{':   '},num2str(round(f(mode),2)),{' [Hz]'});
        title(ax,TITLE,'FontSize',10)
        view(ax,View);      axis(ax,'off');     axis(ax,'equal');
        shading(ax,'interp');  colormap(ax,cmap);  caxis(ax,[-1 1]*max(abs(coorZ)));
    end

end
