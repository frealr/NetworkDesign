close all;
clear all;
clc;

%% Parte 1. Resolvemos el problema para una red SIN Capacitar para el MIP y el MIP regularizado
%Red de 3 nodos
%alpha = 0.5
%utilizo 200 regiones para linealizar la logit en el MIP


cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision high
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 1000.0);
cvx_save_prefs 


num_workers = 8;
% Inicia el parpool (piscina de trabajadores) para parfor
parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar

betas = 1:7; %primero tengo q ver los valores de capacidad que salen para ajustar bien los costes lineales.
budgets = [4e1,6e1,8e1,1e2,1.2e2,1.4e2,1.6e2];


parfor bb=1:length(betas)
    rng(1);
    eps = 1e-3;
    alfa = 1;
    lam = 5;
    beta = betas(bb);

    bud = budgets(bb);
    dm_pax = 1.2e4; 
    dm_op = 1e2; 
    [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim_MIP_entr(lam,bud,alfa,dm_pax,dm_op,beta);
    disp(a);
    disp(f);
    disp(s);
    
    budgets_used_MIP_entr(bb) = budget;
    
    disp(['budget = ',num2str(bud),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
    %disp(f);

end

delete(gcp);

parpool('local', num_workers);



cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision high
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 1000.0);
cvx_save_prefs 

parfor bb=1:length(betas)
    rng(1);
    eps = 1e-3;
    alfa = 1;
    lam = 5;
    beta = betas(bb);
    bud = budgets(bb);
    dm_pax = 1.2e4;
    dm_op = 1e2;

    [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP(lam,bud,alfa,dm_pax,dm_op,beta);
        disp(a);
    disp(f);
    disp(s);

    % disp(['nlinks =',num2str(sum(sum(a > eps)))]);
    % 
    % disp(['budget = ',num2str(budget)]);
     budgets_used_MIP(bb) = budget;
    % 
    % disp(['obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
    %     ', op_obj = ',num2str(op_obj)]);
    % disp(f);
    disp(['budget = ',num2str(bud),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);

end

delete(gcp);

%% See results
lam = 5;
for bb=1:length(betas)
    beta = betas(bb);
    bud = budgets(bb);
    filename = sprintf('./results_paper/sol_MIP_entr_3node_nocap_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
  %  disp(['MIP regularizado, budget = ',num2str(bud),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
    %    ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > 1e-8))),', pres utilizado = ', num2str(budget),', t_comp = ',num2str(comp_time)]);
    obj_val_MIP_entr = obj_val;
    nlinks = sum(sum(a > 1e-8));
    t_comp_entr = comp_time;

    filename = sprintf('./results_paper/sol_MIP_3node_nocap_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
   % disp(['MIP, budget = ',num2str(bud),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
     %   ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > 1e-8))),', pres utilizado = ', num2str(budget),', t_comp = ',num2str(comp_time)]);
    obj_val_MIP = obj_val;
    dif = 100.*(obj_val_MIP_entr-obj_val_MIP)./obj_val_MIP;
   % disp(['dif = ',num2str(dif),' %'])

    disp([num2str(bud),'&',num2str(lam),'&',num2str(obj_val_MIP),'&',num2str(obj_val_MIP_entr),'&',num2str(dif),'&',num2str(budget),'&',num2str(budget),'&', ...
        num2str(nlinks),'&',num2str(nlinks),'&',num2str(comp_time),'&',num2str(t_comp_entr),'\\ \hline']);

end


%% Functions

function budget = get_budget(s,s_prim,a,a_prim,n,...
    station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
    budget = 0;
    for i=1:n
        if s_prim(i) > 1e-8
            budget = budget + lam*station_cost(i); 
        end
        for j=1:n
            if a_prim(i,j) > 1e-8
                budget = budget + lam*link_cost(i,j);
            end
        end
    end
end

function [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op)
    n = 3;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;
    dm_pax = 1.2e4;


    op_obj = op_obj + (sum(sum(op_link_cost.*a_prim))); %operational costs
    for o=1:n
        for d=1:n
            pax_obj = pax_obj + demand(o,d).*fext(o,d);
        end
    end
    obj_val = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
end

function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim_MIP_entr(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_3node_network();

    M = 1e1;
     
    eps = 1e-3;
    tic;
    
    cvx_begin quiet
        variable s(n)
        variable s_prim(n)
        variable s_bin(n) binary
        variable delta_s(n)
        variable a(n,n)
        variable a_prim(n,n)
        variable a_bin(n,n) binary
        variable delta_a(n,n)
        variable f(n,n)
        variable fext(n,n)
        variable fij(n,n,n,n)
        op_obj = 0;
        pax_obj = 0;
        op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operation costs
        for o=1:n
            for d=1:n
                pax_obj = pax_obj + 1e-6*(demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d)))); 
            end
        end
        pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
        pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
        pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
        obj = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
        obj = 1e6.*obj; %quitar para 9 nodos.
        minimize obj
        % constraints
        s >= 0;
        s <= 1e3;
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
        M.*s_bin >= s_prim;
        M.*a_bin >= a_prim;
        bg = 0;
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom.*M; %cuidado, esta sin capacidad
            end
        end
        for i=1:n
            eta*sum(a(:,i)) <= s(i);
            eta*sum(a(i,:)) <= s(i);
        end
        sum(sum(a)) <= a_max;
        for o=1:n
            for d=1:n
                sum(fij(o,:,o,d)) == f(o,d);
            end
        end
        
        for o=1:n
            squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(1 - fext(o,[1:(o-1),(o+1):n])); 
        end
        for d=1:n
            squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -1 + fext([1:(d-1),(d+1):n],d);
        end
        for i=1:n
            fij(i,i,:,:) == 0;
            squeeze(sum(fij(i,:,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),2)) - squeeze(sum(permute(fij(:,i,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),[2,1,3,4]),2)) == 0;
        end
        for o=1:n
            fext(o,o) == 0;
            f(o,o) == 0;
            fij(:,o,o,:) == 0;
        end
        for d=1:n
            fij(d,:,:,d) == 0;
        end
        for o=1:n
            for d=1:n
                if o ~= d
                    f(o,d) + fext(o,d) == 1;
                end
            end
        end
        for i=1:n
            for j=1:n
                if ~ismember(j,candidates{i})
                    a_prim(i,j) == 0;
                end
            end
        end

    cvx_end
    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;

    
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./results_paper/sol_MIP_entr_3node_nocap_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_3node_network();



    M = 1e1;

    nreg = 200;


    eps = 1e-3;
    vals_regs = linspace(0.995,0.005,199);
    [lin_coef,dmax,dmin] = get_linearization(n,nreg,alt_time,alt_price,-1,-1,vals_regs);
    tic;
    
    cvx_begin quiet
        variable s(n)
        variable s_prim(n)
        variable s_bin(n) binary
        variable delta_s(n)
        variable a(n,n)
        variable a_prim(n,n)
        variable a_bin(n,n) binary
        variable delta_a(n,n)
        variable f(n,n)
        variable z(n,n) binary
        variable zij(n,n,n,n)
        variable fext(n,n)
        variable fij(n,n,n,n) binary
        variable delta_rod(nreg,n,n)
        variable alfa_rod(nreg,n,n) binary
        op_obj = 0;
        pax_obj = 0;
        op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operation costs
        for o=1:n
            for d=1:n
                pax_obj = pax_obj + 1e-6*(demand(o,d).*fext(o,d)); 
            end
        end
     %   pax_obj = 1e-6*(demand(1,2).*fext(1,2));
        obj = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
        obj = 1e6.*obj; %quitar para 9 nodos.
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
        zij >= 0;
        zij <= 1;
        zij <= M.*fij;
        s_prim == s + delta_s;
        a_prim == a + delta_a;
        M.*s_bin >= s_prim;
        M.*a_bin >= a_prim;
        bg = 0;
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(zij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom.*M; %cuidado, ahora está sin capacidad
            end
        end


        for i=1:n
            eta*sum(a(:,i)) <= s(i);
            eta*sum(a(i,:)) <= s(i);
        end
        sum(sum(a)) <= a_max;

        
        for o=1:n
            squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(z(o,[1:(o-1),(o+1):n])); 
        end
        for d=1:n
            squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -z([1:(d-1),(d+1):n],d);
        end
        for i=1:n
            fij(i,i,:,:) == 0;
            squeeze(sum(fij(i,:,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),2)) - squeeze(sum(permute(fij(:,i,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),[2,1,3,4]),2)) == 0;
        end
        for o=1:n
            fext(o,o) == 0;
            f(o,o) == 0;
            fij(:,o,o,:) == 0;
        end
        for d=1:n
            fij(d,:,:,d) == 0;
        end
        for o=1:n
            for d=1:n
                if o ~= d
                    f(o,d) + fext(o,d) == 1;
                    f(o,d) <= M.*z(o,d);
                end
            end
        end
        for i=1:n
            for j=1:n
                if ~ismember(j,candidates{i})
                    a_prim(i,j) == 0;
                end
            end
        end

        for o=1:n
            for d=1:n
                zij(:,:,o,d) >= f(o,d) - M.*(1-fij(:,:,o,d));
                zij(:,:,o,d) <= f(o,d).*ones(n);
            end
        end

        delta_rod >= 0;
        for o=1:n
            for d=1:n
                sum(sum(fij(:,:,o,d).*(travel_time + prices))) == sum(delta_rod(:,o,d)); 
            end
        end
        

         for o=1:n
             for d=1:n
                 ss = logit(dmin(1,o,d),-1,-1,alt_time(o,d),alt_price(o,d));
                 for r=2:nreg
                     if dmin(r,o,d) == 0
                         ss = logit(dmin(r,o,d),-1,-1,alt_time(o,d),alt_price(o,d));
                     end
                 end
                 f(o,d) <= ss + sum(delta_rod(:,o,d).*lin_coef(:,o,d));
                 delta_rod(1,o,d) <= dmax(1,o,d) - dmin(1,o,d);
                 for r=2:nreg
                     alfa_rod(r,o,d) <= alfa_rod(r-1,o,d);
                     alfa_rod(r,o,d)*(dmax(r,o,d)-dmin(r,o,d)) <= delta_rod(r,o,d);
                     delta_rod(r,o,d) <= alfa_rod(r-1,o,d)*(dmax(r,o,d)-dmin(r,o,d));
                 end
             end
         end


    cvx_end


    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./results_paper/sol_MIP_3node_nocap_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

end


function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_3node_network()

    n = 3;
    
    %candidates to construct a link for each neighbor
    candidates = {[2,3],[1,3],[1,2]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.5,2;...
                1.5,0,3;...
                2,3,0];
    
    %fixed cost for constructing links
    link_cost = [0,2,5;...
                2,0,1;...
                5,1,0];
    link_cost (link_cost ==0) = 1e4;

    
    %fixed cost for constructing stations
    station_cost = [2,3,2];
    
    link_capacity_slope = link_cost; 
    station_capacity_slope = station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,2,1;...
                   2,0,4;...
                   1,4,0];
    
    distance = [0,1,1;...
                1,0,2;...
                1,2,0];
    
    
    %Load factor on stations
    load_factor = 0.25 .* ones(1, n);
    
    % Op Link Cost
    op_link_cost = 4.*distance;
    
    % Congestion Coefficients
    congestion_coef_stations = 0.1 .* ones(1, n);
    congestion_coef_links = 0.1 .* ones(n);
    
    % Prices
    prices = 0.1.*(distance).^(0.7);
    %prices = zeros(n);
    
    % Travel Time
    travel_time = 60 .* distance ./ 30; % Time in minutes

    
    % Alt Time
    alt_time = 60 .* alt_cost ./ 30; % Time in minutes

    alt_price = 0.1.*(alt_cost).^(0.7); %price

    
    
    a_nom = 588;             
    
    tau = 0.57;
    eta = 0.25;
    a_max = 1e9;
    eps = 1e-3;
end




function [lin_coef,dmax,dmin] = get_linearization(n,nreg,alt_time,alt_price,omega_t,omega_p,vals_regs)
    dmax = zeros(nreg,n,n);
    dmin = dmax;
    lin_coef = dmax;
    
    for o=1:n
        for d=1:n
            syms x
            for r=1:(nreg-1)
                eqn = exp(-x)/( exp(-x) + exp(omega_t*alt_time(o,d) + omega_p*alt_price(o,d)) ) - vals_regs(r) == 0;
                dmax(r,o,d) = solve(eqn,x);
                dmax(r,o,d) = max(0,dmax(r,o,d));
            end
            dmax(nreg,o,d) = 3e1;
            dmin(1,o,d) = 0;
            for r=2:nreg
                dmin(r,o,d) = dmax(r-1,o,d);
            end

            for r=1:nreg
                if (dmax(r,o,d) == dmin(r,o,d))
                    lin_coef(r,o,d) = 0;
                else
                    lin_coef(r,o,d) = (logit(dmax(r,o,d),omega_t,omega_p,alt_time(o,d),alt_price(o,d)) - logit(dmin(r,o,d),omega_t,omega_p,alt_time(o,d),alt_price(o,d))) ./(dmax(r,o,d)-dmin(r,o,d));
                end
            end
        end
    end

end

function val = logit(x,omega_t,omega_p,time,price)
    val = exp(-x)./( exp(-x) + exp(omega_t*time + omega_p*price) );
end