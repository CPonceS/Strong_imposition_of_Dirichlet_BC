function cmap = custom_colormap(m,f)

n = ceil(m/4);
u = [(1:1:n)/n ones(1,n-1) (n:-1:1)/n]';

g = ceil(n/2)-(mod(m,4)==1)+(1:length(u))';  % Green
r = g + n;                                   % Red
b = g - n;                                   % Blue
r(r>m) = []; g(g>m) = [];  b(b<1) = [];      % RGB

cmap = zeros(m,3);
cmap(r,1) = u(1:length(r));                  % Red
cmap(g,2) = u(1:length(g));                  % Green
cmap(b,3) = u(end-length(b)+1:end);          % Blue
cmap = cmap(ceil(m*f):ceil(m*(1-f)),:);

end
