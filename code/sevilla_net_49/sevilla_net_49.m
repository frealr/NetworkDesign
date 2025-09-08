


clear all; close all; clc;

[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,alt_price,a_nom,tau,sigma,...
    a_max,candidates] = parameters_large_sevilla_network();

eps = 1e-3;
alfa = 0.6; %ver que ponemos
dm_pax = 1;
dm_op = 1;

cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision high
cvx_save_prefs 

num_workers = 8;
% Inicia el parpool (piscina de trabajadores) para parfor
%parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
% parpool(num_workers);
lam = 5;
betas = [0,0.01,0.02,0.04,0.06,0.08,0.1,1,10,100,1e3,1e4];

%parfor ll=1:length(lams)
 %   lam = lams(ll);
for bb=1:length(betas)
    eps = 1e-3;
    alfa = 1; %ver que ponemos
    beta = betas(bb);
    dm_pax = 1;
    dm_op = 1;
    [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim(lam,beta,alfa,dm_pax,dm_op);

    disp(['nlinks =',num2str(sum(sum(a_prim > eps)))]);
    
    disp(['budget = ',num2str(budget)]);
    
    disp(['obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj)]);
end

%end


% Cierra el parpool
%delete(gcp);



function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim(lam,beta_or,alfa,dm_pax,dm_op)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,alt_price,a_nom,tau,sigma,...
    a_max,candidates] = parameters_large_sevilla_network();


    niters = 10;           
    eps = 1e-3;
    a_prev = 1e4.*ones(n);
    s_prev= 1e4.*ones(n,1);
    logit_coef = 0.25;
    disp(['beta = ',num2str(beta_or),', lam = ',num2str(lam)]);
    tic;
    for iter=1:niters
        cvx_begin

            variable s(n)
            variable s_prim(n)
            variable delta_s(n)
            variable a(n,n)
            variable a_prim(n,n)
            variable delta_a(n,n)
            variable f(n,n)
            variable fext(n,n)
            variable fij(n,n,n,n)
            op_obj = 0;
            pax_obj = 0;
            bud_obj = 0;
            bud_obj = bud_obj + 1e-3*sum(station_capacity_slope.*s_prim);
            bud_obj = bud_obj + 1e-3*(sum(sum(link_capacity_slope.*a_prim)));  %linear construction costs
            op_obj = op_obj + 1e-3*(sum(sum(op_link_cost.*a_prim))); %operation costs

            if iter < niters
                pax_obj = pax_obj + 1e-3 * sum(inv_pos(congestion_coef_stations .* delta_s + eps));
                pax_obj = pax_obj + 1e-3 * sum(sum(inv_pos(congestion_coef_links .* delta_a + eps)));
                bud_obj = bud_obj + 1e-3*lam*sum(sum((link_cost.*a_prim.*(1./(a_prev+eps))))) + 1e-3*lam*sum((station_cost.*s_prim.*(1./(s_prev+eps)))); %fixed construction costs         
            end
            

            prices_exp = repmat(prices, 1, 1, n, n);  % size [n x n x n x n]
            
            % Compute contribution per (o,d)
            contrib = squeeze(sum(sum(prices_exp .* fij * logit_coef, 1), 2));  % size [n x n]
            
            % Final pax_obj accumulation
            pax_obj = pax_obj + 1e-3 * sum(sum(demand .* contrib));


            pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(alt_price).*fext.*logit_coef)));
            pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(-entr(f) - f))));
            pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(-entr(fext) - fext))));
    
    
            if iter == niters
                for i=1:n
                    if s_prev(i) >= eps
                        pax_obj = pax_obj + 1e-3*(sum(inv_pos(congestion_coef_stations(i).*delta_s(i) + eps)));
                    end
                    for j=1:n
                        if a_prev(i,j) >= eps
                            pax_obj = pax_obj + 1e-3*(sum(sum(inv_pos(congestion_coef_links(i,j).*delta_a(i,j) + eps))));
                         end
                     end
                 end
            end
            obj = beta_or*bud_obj + (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;

            minimize obj
            % constraints
            s >= 0;
            s_prim >= 0;
            delta_s >= 0;
            a >= 0;
            a_prim >= 0;
            delta_a >= 0;
            f >= 0;
            f <= 1;
            fij >= 0;
            fij <= 1;
            fext >= 0;
            fext <= 1;
            s_prim == s + delta_s;
            a_prim == a + delta_a;


            a_prim == a_prim';

         
            % Expand demand to match dimensions of fij: [1 x 1 x n x n]
            demand_expanded = repmat(reshape(demand, [1 1 n n]), [n n 1 1]);  % [n x n x n x n]

            
            % Multiply: broadcasting demand across the first two dimensions
            flow_link = squeeze(sum(sum(fij .* demand_expanded, 3), 4));  % size: [n x n]
            
            % Final constraint matrix:
            flow_link <= tau * a * a_nom;

            sigma * sum(a_prim, 1)' <= s;

            sum(sum(a_prim)) <= a_max;

            S = squeeze(sum(fij,2));
            mask = repmat(eye(n),1,1,n);
            squeeze(sum(S.*mask,2)) == f;

            % máscara para excluir la diagonal
            % fij: n×n×n×n
            % fext: n×n
            E  = eye(n);
            
            % Sumatorios que usabas en el bucle
            S2 = squeeze(sum(fij, 2));   % n×n×n,  S2(i,k,j) = sum_p fij(i,p,k,j)
            S1 = squeeze(sum(fij, 1));   % n×n×n,  S1(o,k,j) = sum_i fij(i,o,k,j)
            
            % Extraer la diagonal k=i (o=k) en cada página j, sin bucles:
            A = squeeze(sum(repmat(E, [1 1 n]) .* S2, 2));   % n×n, A(o,j) = sum_p fij(o,p,o,j)
            B = squeeze(sum(repmat(E, [1 1 n]) .* S1, 2));   % n×n, B(o,j) = sum_i fij(i,o,o,j)
            
            % Residual de las igualdades (A - B) == (1 - fext) solo fuera de la diagonal
            X   = (A - B) - (1 - fext);         % n×n
            Y   = X.';                          % reordenamos para empaquetar por filas o
            res = reshape(Y(~eye(n)), n-1, n)'; % n×(n-1), como en tu bucle: por fila o, j≠o
            
            % Si quieres booleano como en tu '==':
            res == 0;   % o usa tolerancia: abs(res) <= tol



     %       for o=1:n
     %           squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(1 - fext(o,[1:(o-1),(o+1):n])); 
     %       end
          

            % Sumar sobre la 2ª dim y sobre la 1ª dim (como en tu código original)
            S2 = squeeze(sum(fij, 2));   % n×n×n,  S2(i,k,l) = sum_j fij(i,j,k,l)
            S1 = squeeze(sum(fij, 1));   % n×n×n,  S1(j,k,l) = sum_i fij(i,j,k,l)
            
            % Necesitamos, para cada d, los términos con l=d y k≠d:
            % A(d,k) = sum_j fij(d,j,k,d) = S2(d,k,d)
            % B(d,k) = sum_i fij(i,d,k,d) = S1(d,k,d)
            
            n  = size(S2,1);
            [I,K] = ndgrid(1:n, 1:n);            % rejilla de índices (d,k)
            idx = sub2ind([n,n,n], I, K, I);     % (d,k,l=d) en S2/S1
            
            A = S2(idx);                         % n×n con A(d,k) = S2(d,k,d)
            B = S1(idx);                         % n×n con B(d,k) = S1(d,k,d)
            
            % Residuo de la igualdad: (A - B) == (-1 + fext)
            X = (A - B) - (-1 + fext);           % n×n
            
            % Reorganizar a n×(n-1) excluyendo la diagonal k=d (mismo orden que [1:d-1,d+1:n])
            Y   = X.';                            % transponer para aplanar por filas d
            res = reshape(Y(~eye(n)), n-1, n)';   % n×(n-1): por fila d, columnas k≠d
            
            % Si quieres booleano exacto (mejor usar tolerancia):
            % tol = 1e-9;  eqs = abs(res) <= tol;
            res == 0;



        %    for d=1:n
        %        squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -1 + fext([1:(d-1),(d+1):n],d);
        %    end

                    % --- 1) Restricción fij(i,i,:,:) == 0 (diagonal en las dos primeras dims) ---
            E = eye(n);
            mask12 = repmat(E, [1 1 n n]);         % n×n×n×n, true donde i==j
            fij(mask12 == 1) == 0;              % vector lógico con las igualdades de la diagonal
            
            % --- 2) Restricción sum_j fij(i,j,k,l) - sum_i fij(i,j=k,l) == 0 para k≠i, l≠i ---
            % Sumas sobre dims 2 y 1, respectivamente (mismas que en tu for)
            S2 = squeeze(sum(fij, 2));             % n×n×n  con S2(i,k,l) = sum_j fij(i,j,k,l)
            S1 = squeeze(sum(fij, 1));             % n×n×n  con S1(i,k,l) = sum_i fij(i,j,k,l)
            
            res = S2 - S1;                         % n×n×n, residual de la igualdad
            
            % Máscara de posiciones donde imponer la igualdad: k≠i y l≠i
            [ii,kk,ll] = ndgrid(1:n, 1:n, 1:n);
            valid = (kk ~= ii) & (ll ~= ii);       % n×n×n
            
            % Igualdad vectorizada solo en las posiciones válidas
            (res(valid) == 0);               % vector lógico de igualdades válidas
            % (si prefieres un tensor lógico con las no-válidas ignoradas)
            (res .* valid) == 0;      % n×n×n lógico
            %for i=1:n
              %  fij(i,i,:,:) == 0;
             %   squeeze(sum(fij(i,:,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),2)) - squeeze(sum(permute(fij(:,i,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),[2,1,3,4]),2)) == 0;
            %end

               
            fext.*eye(n) == 0;
            f.*eye(n) == 0;


            % fij: n×n×n×n

            n = size(fij,1);
            
            % --- 1) fij(:,o,o,:) == 0  (segunda y tercera dimensión iguales)
            mask1 = false(n,n,n,n);
            for o = 1:n
                mask1(:,o,o,:) = true;
            end
             (fij(mask1) == 0);   % vector lógico con esas igualdades
            
            % --- 2) fij(d,:,:,d) == 0  (primera y cuarta dimensión iguales)
            mask2 = false(n,n,n,n);
            for d = 1:n
                mask2(d,:,:,d) = true;
            end
            (fij(mask2) == 0);
            
            % Si prefieres en un solo tensor lógico (para ver todas las posiciones nulas):
            mask = mask1 | mask2;
            eq_all = (fij(mask) == 0);
            
          


            % for o=1:n
            %     fij(:,o,o,:) == 0;
            % end
            % for d=1:n
            %     fij(d,:,:,d) == 0;
            % end
            f+fext == 1-eye(n);


           


           mask = false(n);                 % n×n lógico
            for i = 1:n
                mask(i, candidates{i}) = true;   % marcamos como válidos
            end
            
            % Los que NO están en candidates{i} para cada fila i:
            forbidden = ~mask;   % matriz n×n con true donde debería valer 0
            
            % Igualdades vectorizadas:
            (a_prim(forbidden) == 0);
            % 
            % for i=1:n
            % 
            %     for j=1:n
            %         if ~ismember(j,candidates{i})
            %             a_prim(i,j) == 0; %ver esto
            %         end
            %     end
            % end
            if iter == niters
                s_prim(s_prev <= 0.1) == 0;
                a_prim(a_prev <= 0.1) == 0;
            end
                       disp('hola')

    
        cvx_end
        a_prev = a_prim;
        s_prev = s_prim;

        disp(['iter = ',num2str(iter),', beta = ',num2str(beta_or),', lam = ',num2str(lam),', nlinks =',num2str(sum(sum(a_prim > 0.1))),', nstations = ', ...
            num2str(sum(s_prim > 0.1)),', att_dem = ',num2str(100*sum(sum(f.*demand))/sum(sum(demand))),' %']);

    end
    comp_time = toc;
    
    
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./results_paper/large_net_prueba_beta=%d_lam=%d.mat',beta_or,lam);

    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

end


function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,alt_price,a_nom,tau,sigma,...
    a_max,candidates] = parameters_large_sevilla_network()


    n_fs = 140; %number of demand centroids
    
    nlines = 8; %resolution for the grid
    sigmaX = 3;
    sigmaY= sigmaX;

    [n,demand,candidates,distance] = set_network(n_fs,...
    nlines,sigmaX,sigmaY,'g');
    
    
    %fixed cost for constructing links

    link_cost = 2e7.*distance./(365.25*25);

    link_capacity_slope = 0.03.*link_cost; 

    station_cost = 1e7./(365.25*25);

    station_capacity_slope = 0.02.*station_cost;
    
    
    % Op Link Cost
    op_link_cost = 3.*distance;
    
    % Congestion Coefficients
    congestion_coef_stations = 0.1 .* ones(n, 1);
    congestion_coef_links = 0.1 .* ones(n);
    
    % Prices
    riding_cost = 0.083; %e/min
    train_speed = 40; %km/h
    fare = 0; %e

    

    prices = fare + (riding_cost*train_speed/60).*distance;

    %Alternative: car
    
    fixed_alt_price = 1.75; %e
    variable_alt_price = 0.12; %e/km
    alt_time_val = 0.05; %e/min
    
    average_parking_time = 10; %min
    parking_cost = 0.25; %e/min
    fixed_parking_cost = 1; %e
    alt_speed = 60; %km/h
    
    alt_price = fixed_alt_price + variable_alt_price.*distance ...
    + alt_time_val*(alt_speed/60).*distance ...
    + fixed_parking_cost + parking_cost*average_parking_time;
    
    a_nom = 607;             
    
    tau = 0.57;
    sigma = 0.25;
    a_max = 1e9;
    eps = 1e-3;
end




function [n,demand_cells_sub,candidates,distance] = set_network(n_fs,...
    nlines,sigmaX,sigmaY,cell_distribution)

    coordinates = readtable('./NODOSSEV_MT.xlsx','Sheet',1);
    id_node = table2array(coordinates(1:140,1));
    coor_x_raw = table2array(coordinates(1:140,2));
    coor_y_raw = table2array(coordinates(1:140,3));
    
    %real distance between nodes 1 and 60
    lat1 = 37.413969; lon1= -6.007699;
    lat60 = 37.361042; lon60 = -5.959065;
    real_d = haversine(lat1,lon1,lat60,lon60);
    
    %scaled distance between nodes 1 and 60
    scaled_d = sqrt( (coor_x_raw(id_node == 1) - coor_x_raw(id_node == 60))^2 ...
        +(coor_y_raw(id_node == 1) - coor_y_raw(id_node == 60) )^2   );
    
    scale = real_d / scaled_d; % scale of coordinates
    
    coor_x(id_node) = coor_x_raw*scale;
    coor_y(id_node) = coor_y_raw*scale;
    
    % constant values
    
    meters_per_deg_lat = 111.320; 
    meters_per_deg_lon = 111.320 * cosd(lat1);  
    
    % distance wrt reference node
    dx = coor_x - coor_x(id_node==1);
    dy = coor_y - coor_y(id_node==1);
    
    % Approximate lat/lon
    lat = lat1 + dy / meters_per_deg_lat;
    lon = lon1 + dx / meters_per_deg_lon;
    lon(18) = lon(18) + 0.0005;
    lat(18) = lat(18) + 0.0005;
    
    lon(34) = lon(34) + 0.0005;
    lat(34) = lat(34) + 0.0005;
    %
    
    demand_raw = readtable('./NODOSSEV_MT.xlsx','Sheet',2);
    id_node_x = table2array(demand_raw(1,2:141));
    id_node_y = table2array(demand_raw(2:141,1));
    
    demand_aux = table2array(demand_raw(2:141,2:141));
    
    demand(id_node_x,id_node_y) = demand_aux;
    
    %full adyacency matrix
    A_full = zeros(n_fs);
    G_full = graph(A_full);
    figure;
    
    scatter(coor_x,coor_y,'filled');
    hold on;
    xlims = xlim; ylims = ylim;
    
    ix = (1:nlines) - 0.5 - nlines/2;
    iy = (1:nlines) - 0.5 - nlines/2;
    
    densX = exp(-(ix.^2)/(2*sigmaX^2));
    densY = exp(-(iy.^2)/(2*sigmaY^2));
    
    weights = @(n,sigma)...
        exp(-(( (1:n) -0.5 - n/2  ).^2 ) / (2*sigma^2));
    normalize = @(v) v / sum(v);

    switch cell_distribution
        case 'u'
            % uniform split for cells
            wx = ones(1,nlines)/nlines;
            wy = ones(1,nlines)/nlines;
        otherwise
            % gaussian split for cells
            wx = normalize(1./densX);
            wy = normalize(1./densY);
    end
    
    widths = wx * (xlims(2) - xlims(1));
    heights = wy * (ylims(2) - ylims(1));
    
    xedges = [xlims(1),xlims(1) + cumsum(widths)];
    
    yedges =  [ylims(1),ylims(1) + cumsum(heights)];
    
    x_cent = (xedges(1:end-1)+xedges(2:end))/2;
    y_cent = (yedges(1:end-1)+yedges(2:end))/2;
    
    [Xc,Yc] = meshgrid(x_cent,y_cent);
    
    xc = Xc(:);
    yc = Yc(:);
    
    for xv=xedges
        xline(xv,':k','LineWidth',0.8);
    end
    for yv=yedges
        yline(yv,':k','LineWidth',0.8);
    end
    
    
    N = numel(coor_x);
    
    ix = discretize(coor_x, xedges);  
    iy = discretize(coor_y, yedges); 
    
    valid = ~isnan(ix) & ~isnan(iy);     
    cell_lin = nan(N,1);
    cell_lin(valid) = sub2ind([nlines nlines], iy(valid), ix(valid));   
    
    K = nlines*nlines;  
    
    
    demanda_total = sum(demand, 2) + sum(demand, 1)';
    
    centroides_x = nan(K,1);
    centroides_y = nan(K,1);
    
    for k = 1:K
        % nodes in cell k
        idx_nodos = find(cell_lin == k);
    
        if isempty(idx_nodos)
            centroides_x(k) = Xc(k);
            centroides_y(k) = Yc(k);
            continue;  % empty cell
        end
    
        % total demand per cell
        dt = demanda_total(idx_nodos);
    
        % largest demand node index per cell, allocate station there
        [~, idx_max] = max(dt);
        idx_nodo_max = idx_nodos(idx_max);
    
        centroides_x(k) = lat(idx_nodo_max);
        centroides_y(k) = lon(idx_nodo_max);
    end
    
    [I,J,V] = find(demand);                   
    ci = cell_lin(I);  cj = cell_lin(J);  
    
    mask = ~isnan(ci) & ~isnan(cj);      
    demand_cells = accumarray([ci(mask), cj(mask)], V(mask), [K K], @sum, 0);
    demand_cells(1:size(demand_cells,1)+1:end) = 0;
    
    
    nodos_inactivos = 0;
    for i=1:length(demand_cells)
        if sum(demand_cells(i,:)) == 0
            nodos_inactivos = nodos_inactivos + 1;
        end
    end
    
    n = size(Xc, 1);            
    K = n * n;
    
    % Initialize candidate links matrix
    A = zeros(K, K);
    
    % Connection probability
    p = 0.8;
    rng(123);
    
    for row = 1:n
        for col = 1:n
            i = sub2ind([n n], row, col); 
    
            neighbors = [];
    
            if row > 1     
                j = sub2ind([n n], row-1, col);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1; 
                end
            end
            if row < n    
                j = sub2ind([n n], row+1, col);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1;
                end
            end
            if col > 1  
                j = sub2ind([n n], row, col-1);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1;
                end
            end
            if col < n 
                j = sub2ind([n n], row, col+1);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1;
                end
            end
        end
    end
    
    xpos = centroides_x;   
    ypos = centroides_y;   
    G = graph(A);
    
    % aggregate demand per cell
    demanda_por_celda = accumarray(cell_lin(valid), demanda_total(valid), [K 1], @sum, 0);
    dtotal = demand_cells + demand_cells';
    zonas_con_demanda = find(sum(dtotal) > 0);
    G_sub = subgraph(G, zonas_con_demanda);
    
    % potential link matrix
    A_sub = full(adjacency(G_sub));
    demand_cells_sub = demand_cells(zonas_con_demanda,zonas_con_demanda);
    
    % coordinates of cells with demand
    x_sub = xpos(zonas_con_demanda);
    y_sub = ypos(zonas_con_demanda);
    
    figure; hold off; axis equal;
    geoscatter(lat, lon, 40, 'r', 'filled');  
    hold on;
    
    
    [i_idx, j_idx] = find(A_sub); 
    linecolor = [0, 0, 0, 0.2];
    for k = 1:length(i_idx)
        i = i_idx(k);
        j = j_idx(k);
        geoplot([x_sub(i), x_sub(j)], [y_sub(i), y_sub(j)], 'k-', 'LineWidth', 1.5,'Color',linecolor); hold on;
    end
    
    for n = 1:numel(x_sub)
        text(x_sub(n), y_sub(n), sprintf('%d', n), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', ...
             'FontSize', 11, ...
             'Color', 'black');
    end
    
    geobasemap('topographic');
    n = length(x_sub);
    
    coor_x_sub = centroides_x(zonas_con_demanda);
    coor_y_sub = centroides_y(zonas_con_demanda);
    
    K = size(A_sub, 1);
    candidates = cell(K, 1);
    
    for i = 1:K
        candidates{i} = find(A_sub(i, :) == 1);
    end
    distance = 1e6.*ones(n);
    for i=1:n
        distance(i,i) = 0;
        cand = candidates(i);
        cand = cand{1};
        cand = cand(cand > i);
        for j=i+1:n
            if sum(j == cand) > 0
                distance(i,j) = haversine(coor_x_sub(i), coor_y_sub(i), coor_x_sub(j), coor_y_sub(j));
                distance(j,i) = distance(i,j);
            end
        end
    end

end

function distancia = haversine(lat1, lon1, lat2, lon2)
    % Convierte las coordenadas de grados a radianes
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Diferencias en coordenadas
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;

    % Fórmula haversine
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));

    % Radio de la Tierra en kilómetros (aproximado)
    radio_tierra = 6371;

    % Calcula la distancia
    distancia = radio_tierra * c;
end