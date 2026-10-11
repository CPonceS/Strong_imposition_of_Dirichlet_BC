function modesSpatialTimBeam(model,method,N_modes,View)
%
% Plots the first N_modes eigenmodes of a Spatial Timoshenko Beam
%
%   model   : output of femPHS1D
%   method  : 'sFEM' or 'mFEM'
%   N_modes : number of modes to be plotted
%   View    : [Az El] (optional, default [35 30])
%
% State vector (its dimension depends on the model):
%   sFEM:  x = [ p (N_free) | r_free (N_free) | r_D (N_Dir) ]
%   mFEM:  x = [ p (N_free) | eps (M_Omega) | r_free (N_free) | r_D (N_Dir) ]
%
% The displacements are the last N_Omega entries, ordered as [r_free; r_D].
% model.idx_sort puts them in nodal order.
% The last 2 PDE variables are the displacements in the y and z axes.

    if nargin < 4 || isempty(View); View = [35 30]; end

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
    
    % Limits
    lim_x  = [min(coor) max(coor)];
    lim_yz = [-w_max w_max]*1.1;          
    
    wall_y = ones(size(coor))*lim_yz(2);   % XZ
    floor_z = ones(size(coor))*lim_yz(1);  % XY
    
    % Single figure: square or rectangle with the fewest rows that fits N_modes
    nc = ceil(sqrt(N_modes));                % columns
    nr = ceil(N_modes/nc);                   % rows   
    figure('Color','w','Units','normalized','Position',[0.15 0.1 0.7 0.75]);
    
    for mode = 1:N_modes
        r = V(idx_r,idx(mode));
        r = r(model.idx_sort);                         % nodal ordering
        
        % Extract transverse displacements in y and z axes
        wy = r((n_PDE-2)*NN+1:(n_PDE-1)*NN);
        wz = r((n_PDE-1)*NN+1:n_PDE*NN);
        
        % Find the global maximum to properly align the phase for BOTH directions
        [~, imax_y] = max(abs(wy)); max_y = abs(wy(imax_y));
        [~, imax_z] = max(abs(wz)); max_z = abs(wz(imax_z));
        
        if max_y >= max_z
            phase = angle(wy(imax_y));
        else
            phase = angle(wz(imax_z));
        end
        
        % Align the phase -> real mode shapes
        wy = real(wy*exp(-1i*phase));
        wz = real(wz*exp(-1i*phase));
        
        % Normalize maintaining relative magnitude between y and z
        w_norm = max(sqrt(wy.^2 + wz.^2));
        wy_mod = (wy/w_norm)*w_max;
        wz_mod = (wz/w_norm)*w_max;
        
        ax = subplot(nr,nc,mode);                      % row-wise: left to right, top to bottom
        TITLE = strcat({'Mode '},num2str(mode),{':   '},num2str(round(f(mode),2)),{' [Hz]'});
        
        hold(ax,'on'); box(ax,'on');
        
        % XY shadows
        plot3(ax, coor, wy_mod, floor_z, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5)
        plot3(ax, coor, 0*wy_mod, floor_z, 'Color', [0.8 0.8 0.8], 'LineStyle', '--')
        
        % XZ shadows
        plot3(ax, coor, wall_y, wz_mod, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5)
        plot3(ax, coor, wall_y, 0*wz_mod, 'Color', [0.8 0.8 0.8], 'LineStyle', '--')
        
        % 3. Curva principal y línea de referencia
        plot3(ax, coor, wy_mod, wz_mod, 'k', 'LineWidth', 2) 
        plot3(ax, coor, 0*wy_mod, 0*wz_mod, 'k--')           
        
        hold(ax,'off');
        
        title(ax, TITLE, 'FontSize', 10);   
        axis(ax,'equal'); 
        grid(ax,'on');
        view(ax,View);       
        
        xlim(ax,lim_x);
        ylim(ax,lim_yz);
        zlim(ax,lim_yz);
        
        xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    end
end