function modesPlanarTimBeam(model,method,N_modes)
%
% Plots the first N_modes eigenmodes of a Planar Timoshenko Beam
%
%   model   : output of femPHS1D
%   method  : 'sFEM' or 'mFEM'
%   N_modes : number of modes to be plotted
%
% State vector (its dimension depends on the model):
%   sFEM:  x = [ p (N_free) | r_free (N_free) | r_D (N_Dir) ]
%   mFEM:  x = [ p (N_free) | eps (M_Omega) | r_free (N_free) | r_D (N_Dir) ]
%
% The displacements are the last N_Omega entries, ordered as [r_free; r_D].
% model.idx_sort puts them in nodal order: r = [psi (NN); w (NN); ...].

    A = full(model.(method).J*model.(method).Q);

    NN      = model.NN;
    N_free  = model.N_free;
    N_Omega = model.N_Omega;
    M_Omega = model.M_Omega;
    n_PDE   = N_Omega/NN;
    coor    = model.coor;

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

    w_max = abs(max(coor)-min(coor))/4;      % plotted amplitude

    % Single figure: square or rectangle with the fewest rows that fits N_modes
    nc = ceil(sqrt(N_modes));                % columns
    nr = ceil(N_modes/nc);                   % rows   
    figure('Color','w','Units','normalized','Position',[0.15 0.1 0.7 0.75]);

    for mode = 1:N_modes

        r = V(idx_r,idx(mode));
        r = r(model.idx_sort);                         % nodal ordering
        w = r((n_PDE-1)*NN+1:n_PDE*NN);                % transverse displacement

        [~,imax] = max(abs(w));
        w = real(w*exp(-1i*angle(w(imax))));           % align the phase -> real mode shape
        w_mod = (w/max(abs(w)))*w_max;

        subplot(nr,nc,mode)                            % row-wise: left to right, top to bottom
        TITLE = strcat({'Mode '},num2str(mode),{':   '},num2str(round(f(mode),2)),{' [Hz]'});
        hold on , box on
        plot(coor,w_mod,'k','LineWidth',2)
        plot(coor,0*w_mod,'k--')
        hold off
        title(TITLE,'FontSize',10);   axis equal; set(gca,'ytick',[0])
    end

end