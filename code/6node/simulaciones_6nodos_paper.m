close all;
clear all;
clc;

%% Resolvemos los problemas para distintos betas y lams, tenemos el convexo,
% el MIP con la entropía y el MIP con la logit en restricciones linealizado

%Problema del operador:

%op_obj = operation_costs
%obj_val = 0


%Problema de los pasajeros:
%pax_obj = (distances + prices) + entropies + congestion
%obj_val. Se calcula para un presupuesto dado. En este caso daremos un
%valor a beta, calculamos el presupuest obtenido y lo mantenemos todo el
%tiempo.

betas = [0.1,0.25,0.5,1,1.5,2,2.5]; %estos son los de verdad
betas = [2.75,3,3.25,3.5,3.75,4,4.25,4.5];
betas = [5,7.5,10,12.5,15];
betas = [17.5,20,22.5,25];

betas = [7.5,10,12.5,15,17.5,20];

%betas = [2.75,3];
%betas = [2.75,3,3.25,3.5,3.75,4,4.25,4.5];
%betas = [0.1,0.25,1];
budgets = zeros(length(betas),1);
 num_workers = 8;
%%

cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision high
cvx_save_prefs 
% 

% Inicia el parpool para parfor
% parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
% 
% 
% parfor bb=1:length(betas)
%     eps = 1e-3;
%     alfa = 0.5;
%     lam = 5;
%     beta = betas(bb);
%     dm_pax = 0.35;
%     dm_op = 0.008;
%     [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
%     pax_obj,op_obj] = compute_sim(lam,beta,alfa,dm_pax,dm_op);
% 
% 
% 
%     disp(['nlinks =',num2str(sum(sum(a_prim > eps)))]);
% 
%  %   disp(['budget = ',num2str(budget)]);
%     budgets(bb) = budget;
% 
%   %  disp(['obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
%    %     ', op_obj = ',num2str(op_obj)]);
% 
% end


% Cierra el parpool
 % delete(gcp);

% %% Simulo MIPs con entropía con ese presupuesto
% %tener cuidado, hay que esperar a que acaben los primeros para poder lanzar
% %los segundos.
%% 
cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision default
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 3600.0);
cvx_save_prefs 

% Inicia el parpool (piscina de trabajadores) para parfor
%parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
betas = 1:8;
%budgets = 1e4.*[3.4,3,2.5,2,1.5,1];
budgets = 1e4.*linspace(1.5,4,8);

for bb=2:length(betas)
    eps = 1e-3;
    alfa = 0.5;
    lam = 5;
    beta = betas(bb);
    bud = budgets(bb);
    dm_pax = 0.01;
    dm_op = 0.008;
    [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim_MIP_entr(lam,bud,alfa,dm_pax,dm_op,beta);

    budgets_used_MIP_entr(bb) = budget;
    
    disp(['bb = ',num2str(bb),', budget = ',num2str(bud),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);

end
%delete(gcp);
%%
clc;
lam = 5;
for bb=5:length(betas)
    beta = betas(bb);
    filename  = sprintf('./results_paper/sol_prueba_MIP_entr_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    disp(['MIP regularizado: beta = ',num2str(beta), ', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
end

%%


% %% Simulo MIPs sin entropía con ese presupuesto
% %tener cuidado, hay que esperar a que acaben los primeros para poder lanzar
% %los segundos.
% 
cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision default
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 1800.0);
cvx_save_prefs 


% Inicia el parpool (piscina de trabajadores) para parfor
% parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
% 
% for bb=1:length(betas)
%     eps = 1e-3;
%     alfa = 0.5;
%     lam = 5;
%     beta = betas(bb);
%     bud = budgets(bb);
%     dm_pax = 1; 
%     dm_op = 0.008;
%     [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
%     pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_30(lam,bud,alfa,dm_pax,dm_op,beta);
% 
%     budgets_used_MIP_30min(bb) = budget;
% 
%     disp(['budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
%         ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
% 
% 
% end
% 
% delete(gcp);


%%


% %% Simulo MIPs sin entropía con ese presupuesto
% %tener cuidado, hay que esperar a que acaben los primeros para poder lanzar
% %los segundos.
% 
% cvx_solver_settings -clearall
% cvx_solver mosek
% cvx_precision default
% cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 3600.0);
% cvx_save_prefs 
% 
% 
% 
% % Inicia el parpool (piscina de trabajadores) para parfor
% % parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
% % 
% for bb=1:length(betas)
%     eps = 1e-3;
%     alfa = 0.5;
%     lam = 5;
%     beta = betas(bb);
%     bud = budgets(bb);
%     dm_pax = 1; 
%     dm_op = 0.008;
%     [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
%     pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_60(lam,bud,alfa,dm_pax,dm_op,beta);
% 
%     budgets_used_MIP_1h(bb) = budget;
% 
%     disp(['budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
%         ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
% 
% 
% end
% 
% delete(gcp);

%%


% %% Simulo MIPs sin entropía con ese presupuesto
% %tener cuidado, hay que esperar a que acaben los primeros para poder lanzar
% %los segundos.
% 
cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision default
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 28800.0);
cvx_save_prefs 


% Inicia el parpool (piscina de trabajadores) para parfor
% parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
% 
% for bb=1:length(betas)
%     eps = 1e-3;
%     alfa = 0.5;
%     lam = 5;
%     beta = betas(bb);
%     bud = budgets(bb);
%     dm_pax = 1; 
%     dm_op = 0.008;
%     [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
%     pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_240(lam,bud,alfa,dm_pax,dm_op,beta);
% 
%     budgets_used_MIP_4h(bb) = budget;
% 
%     disp(['budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
%         ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
% 
% 
% end

% delete(gcp);


%%


% %% Simulo MIPs sin entropía con ese presupuesto
% %tener cuidado, hay que esperar a que acaben los primeros para poder lanzar
% %los segundos.
% 
cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision default
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 86400.0);
cvx_save_prefs 


% Inicia el parpool (piscina de trabajadores) para parfor
% parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
betas=1:8;
for bb=5:length(betas)
    eps = 1e-3;
    alfa = 0.5;
    lam = 5;
    beta = betas(bb);
    bud = budgets(bb);
    dm_pax = 1; 
    dm_op = 0.008;
    [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_24h(lam,bud,alfa,dm_pax,dm_op,beta);

    budgets_used_MIP_24h(bb) = budget;

    disp(['budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);


end

%% See results
%clc;
lam = 5;
betas = [0.1,0.25,0.5,1,1.5,2,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,5,7.5,10,12.5,15,17.5,20,22.5,25]; %estos son los de verdad
%betas = [7.5,10];
%betas = 1:6;
betas = 1:8;

for bb=1:length(betas)
    beta = betas(bb);
    filename  = sprintf('./results_paper/sol_prueba_MIP_entr_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    disp(['MIP regularizado: beta = ',num2str(beta), ', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
    obj_entr = obj_val;
    
    filename  = sprintf('./results_paper/sol_prueba_MIP_30min_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    disp(['MIP 30 min : beta = ',num2str(beta), ', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', opt_gap = ',num2str(opt_gap), ...
        ' %, pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
    obj_MIP = obj_val;
    disp(['beta = ', num2str(beta),', diferencia 30 min = ',num2str(100*(obj_entr-obj_MIP)./obj_MIP),' %']);
   filename  = sprintf('./results_paper/sol_prueba_MIP_1h_6node_beta=%d_lam=%d.mat',beta,lam);
   load(filename);
   disp(['MIP 1h min : beta = ',num2str(beta), ', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', opt_gap = ',num2str(opt_gap), ...
       ' %, pax_obj = ',num2str(pax_obj), ...
       ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
   obj_MIP = obj_val;
   disp(['beta = ', num2str(beta),', diferencia 1h min = ',num2str(100*(obj_entr-obj_MIP)./obj_MIP),' %']);

    filename  = sprintf('./results_paper/sol_prueba_MIP_4h_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    disp(['MIP 8h : beta = ',num2str(beta), ', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', opt_gap = ',num2str(opt_gap), ...
        ' %, pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
    obj_MIP = obj_val;
    disp(['beta = ', num2str(beta),', diferencia 8h = ',num2str(100*(obj_entr-obj_MIP)./obj_MIP),' %']);

    filename  = sprintf('./results_paper/sol_prueba_MIP_24h_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    disp(['MIP 24h : beta = ',num2str(beta), ', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', opt_gap = ',num2str(opt_gap), ...
        ' %, pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
    obj_MIP = obj_val;
    disp(['beta = ', num2str(beta),', diferencia 24h = ',num2str(100*(obj_entr-obj_MIP)./obj_MIP),' %']);
end

%%
clc;
for bb=1:length(betas)
    beta = betas(bb);

    filename  = sprintf('./results_paper/sol_prueba_MIP_entr_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    avail_bud = budget;

    obj_val_MIP_entr = obj_val;
    opt_gap_entr = opt_gap;
    bud_entr = budget;
    nlinks_entr = sum(sum(a > 1e-8));
    t_comp_entr = comp_time;


    filename  = sprintf('./results_paper/sol_prueba_MIP_30min_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    obj_val_MIP_30 = obj_val;
    opt_gap_MIP_30 = opt_gap;
    bud_MIP_30 = budget;
    nlinks_30 = sum(sum(a > 1e-8));
    t_comp_MIP_30 = comp_time;

    filename  = sprintf('./results_paper/sol_prueba_MIP_1h_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    obj_val_MIP_1h = obj_val;
    opt_gap_MIP_1h = opt_gap;
    bud_MIP_1h = budget;
    nlinks_1h = sum(sum(a > 1e-8));
    t_comp_MIP_1h = comp_time;

    filename  = sprintf('./results_paper/sol_prueba_MIP_4h_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    obj_val_MIP_8h = obj_val;
    opt_gap_MIP_8h = opt_gap;
    bud_MIP_8h = budget;
    nlinks_8h = sum(sum(a > 1e-8));
    t_comp_MIP_8h = comp_time;

    filename  = sprintf('./results_paper/sol_prueba_MIP_24h_6node_beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    obj_val_MIP_24h = obj_val;
    opt_gap_MIP_24h = opt_gap;
    bud_MIP_24h = budget;
    nlinks_24h = sum(sum(a > 1e-8));
    t_comp_MIP_24h = comp_time;


    dif_30 = 100.*(obj_val_MIP_entr-obj_val_MIP_30)./obj_val_MIP_30;
    dif_1h = 100.*(obj_val_MIP_entr-obj_val_MIP_1h)./obj_val_MIP_1h;
    dif_8h = 100.*(obj_val_MIP_entr-obj_val_MIP_8h)./obj_val_MIP_8h;
    dif_24h = 100.*(obj_val_MIP_entr-obj_val_MIP_24h)./obj_val_MIP_24h;
   % disp(['dif = ',num2str(dif),' %'])

    disp([sprintf('%.1e', avail_bud),'&',num2str(lam),'&',sprintf('%.2f',dif_30),'&',...
        sprintf('%.2f',dif_1h),...%'&',sprintf('%.2f',dif_8h),
        '&',sprintf('%.2f',dif_24h),'&',...
        sprintf('%.2f',-opt_gap_MIP_30), ...
        '&',sprintf('%.2f',-opt_gap_MIP_1h),'&', ...%sprintf('%.2f',-opt_gap_MIP_8h),'&',...
        sprintf('%.2f',-opt_gap_MIP_24h),'&','0','&', ...
        sprintf('%.2e',t_comp_MIP_30),'&',sprintf('%.2e',t_comp_MIP_1h),'&',...sprintf('%.2e',t_comp_MIP_8h),...'&',
        sprintf('%.2e',t_comp_MIP_24h),'&',sprintf('%.2e',t_comp_entr),'\\ \hline']);

end

%%
%see values of alt_time
[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_6node_network();

disp( exp(-travel_time(2,1)  - prices(2,1)) ./(exp(-alt_time(2,1)-alt_price(2,1)) + exp(-travel_time(2,1)  - prices(2,1))));


%% Functions

function budget = get_budget(s,s_prim,a,a_prim,n,...
    station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
    budget = 0;
    for i=1:n
        if s_prim(i) > 1e-3
            budget = budget + lam*station_cost(i) + ...
                station_capacity_slope(i)*s_prim(i);
        end
        for j=1:n
            if a_prim(i,j) > 1e-3
                budget = budget + lam*link_cost(i,j) + ...
                    link_capacity_slope(i,j) * a_prim(i,j);
            end
        end
    end
end

function [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op)
    n = 6;

    dm_pax = 1;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;


    op_obj = op_obj + (sum(sum(op_link_cost.*a_prim))); %operational costs
    for o=1:n
        for d=1:n
            pax_obj = pax_obj + demand(o,d).*fext(o,d);
        end
    end
    obj_val = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim(lam,beta_or,alfa,dm_pax,dm_op)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_6node_network();

    niters = 15;           
    M = 1e6;
    eps = 1e-3;
    a_prev = 1e4.*ones(n);
    s_prev= 1e4.*ones(n,1);
    disp(['beta = ',num2str(beta_or),', lam = ',num2str(lam)]);
    tic;
    for iter=1:niters
      
        cvx_begin quiet
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
            bud_obj = bud_obj + 1e-6*sum(station_capacity_slope'.*s_prim);
            bud_obj = bud_obj + 1e-6*(sum(sum(link_capacity_slope.*a_prim)));  %linear construction costs
            op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operation costs
            if iter < niters
                bud_obj = bud_obj + 1e-6*lam*sum(sum((link_cost.*a_prim.*(1./(a_prev+eps))))) + 1e-6*lam*sum((station_cost'.*s_prim.*(1./(s_prev+eps)))); %fixed construction costs  
            end
    
            for o=1:n
                for d=1:n
                    pax_obj = pax_obj + 1e-6*(demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d)))); 
                end
            end
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
    
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
            s_prim == s+ delta_s;
            a_prim == a+ delta_a;
    
            for i=1:n
                for j=1:n
                    tau.*squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom; %cuidado, esta sin capacidad
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
                squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(f(o,[1:(o-1),(o+1):n])); 
            end
            for d=1:n
                squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -f([1:(d-1),(d+1):n],d);
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
            if iter == niters
                for i=1:n
                    if s_prev(i) <= 0.1
                        s_prim(i) == 0;
                    end
                    for j=1:n
                        if a_prev(i,j) <= 0.1
                            a_prim(i,j) == 0;
                         end
                     end
                 end
             end
    
        cvx_end
        a_prev = a_prim;
        s_prev = s_prim;
        
    end
    comp_time = toc;
    
    
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./results_paper/sol_prueba_beta=%d_lam=%d.mat',beta_or,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

end



function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim_MIP_entr(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_6node_network();

    M = 1e4;
     
    eps = 1e-3;
    tic;
    
    cvx_begin 
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
        M.*s_bin >= s_prim;
        M.*a_bin >= a_prim;
        bg = 0;
        bg = bg + sum(station_capacity_slope'.*s_prim);
        bg = bg + sum(sum(link_capacity_slope.*a_prim));
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom; %cuidado, esta sin capacidad
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
    opt_gap = 100.*(1e6.*cvx_optbnd-obj_val)./obj_val;
    filename = sprintf('./results_paper/sol_prueba_MIP_entr_6node_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget','opt_gap');

end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_30(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_6node_network();

    %esto luego hay que quitarlo, lo de alfa = 1
    %alfa = 1;

    M = 1e4;
    nreg = 10;
    %nreg = 7;

    eps = 1e-3;
    %vals_regs = [0.75,0.5,0.25];
    vals_regs = [0.95,0.9,0.87,0.85,0.83,0.8,0.77,0.75,0.73,0.7,0.67,0.65,0.6,0.55,0.5,0.4,0.3,0.2,0.1];
    
    %vals_regs = linspace(0.95,0.05,6);
    vals_regs = linspace(0.97,0.02,nreg-1);
    [lin_coef,dmax,dmin] = get_linearization(n,nreg,alt_time,alt_price,-1,-1,vals_regs);
    tic;
    
    cvx_begin 
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
        bg = bg + sum(station_capacity_slope'.*s_prim);
        bg = bg + sum(sum(link_capacity_slope.*a_prim));
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(zij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom; %cuidado, ahora está sin capacidad
            end
        end


        for i=1:n
            eta*sum(a(:,i)) <= s(i);
            eta*sum(a(i,:)) <= s(i);
        end
        sum(sum(a)) <= a_max;
        % for o=1:n
        %     for d=1:n
        %         sum(fij(o,:,o,d)) == z(o,d);
        %     end
        % end
        
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

    %tengo que poner las restricciones de la logit, sacar una funcion para
    %obtener las pendientes en funcion de omega y text. probablemente el dm
    %cambie para el problema linealizado, estudiar si lo obtengo o resuelvo
    %para alfa = 1.
    cvx_end


    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;

    %disp(alfa_rod(:,1,2)); disp(delta_rod(:,1,2));
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    opt_gap = 100.*(1e6.*cvx_optbnd-obj_val)./obj_val;
    %disp(f); disp(a);
    filename = sprintf('./results_paper/sol_prueba_MIP_30min_6node_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget','opt_gap');

end



function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_60(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_6node_network();

    %esto luego hay que quitarlo, lo de alfa = 1
    %alfa = 1;

    M = 1e4;
    nreg = 10;
    %nreg = 7;

    eps = 1e-3;
    %vals_regs = [0.75,0.5,0.25];
    vals_regs = [0.95,0.9,0.87,0.85,0.83,0.8,0.77,0.75,0.73,0.7,0.67,0.65,0.6,0.55,0.5,0.4,0.3,0.2,0.1];
    
    %vals_regs = linspace(0.95,0.05,6);
    vals_regs = linspace(0.97,0.02,nreg-1);
    [lin_coef,dmax,dmin] = get_linearization(n,nreg,alt_time,alt_price,-1,-1,vals_regs);
    tic;
    
    cvx_begin 
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
        bg = bg + sum(station_capacity_slope'.*s_prim);
        bg = bg + sum(sum(link_capacity_slope.*a_prim));
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(zij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom; %cuidado, ahora está sin capacidad
            end
        end


        for i=1:n
            eta*sum(a(:,i)) <= s(i);
            eta*sum(a(i,:)) <= s(i);
        end
        sum(sum(a)) <= a_max;
        % for o=1:n
        %     for d=1:n
        %         sum(fij(o,:,o,d)) == z(o,d);
        %     end
        % end
        
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

    %tengo que poner las restricciones de la logit, sacar una funcion para
    %obtener las pendientes en funcion de omega y text. probablemente el dm
    %cambie para el problema linealizado, estudiar si lo obtengo o resuelvo
    %para alfa = 1.
    cvx_end


    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;

    %disp(alfa_rod(:,1,2)); disp(delta_rod(:,1,2));
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    opt_gap = 100.*(1e6.*cvx_optbnd-obj_val)./obj_val;
    %disp(f); disp(a);
    filename = sprintf('./results_paper/sol_prueba_MIP_1h_6node_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget','opt_gap');

end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_240(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_6node_network();

    %esto luego hay que quitarlo, lo de alfa = 1
    %alfa = 1;

    M = 1e4;
    nreg = 10;
    %nreg = 7;

    eps = 1e-3;
    %vals_regs = [0.75,0.5,0.25];
    vals_regs = [0.95,0.9,0.87,0.85,0.83,0.8,0.77,0.75,0.73,0.7,0.67,0.65,0.6,0.55,0.5,0.4,0.3,0.2,0.1];
    
    %vals_regs = linspace(0.95,0.05,6);
    vals_regs = linspace(0.97,0.02,nreg-1);
    [lin_coef,dmax,dmin] = get_linearization(n,nreg,alt_time,alt_price,-1,-1,vals_regs);
    tic;
    
    cvx_begin 
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
        bg = bg + sum(station_capacity_slope'.*s_prim);
        bg = bg + sum(sum(link_capacity_slope.*a_prim));
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(zij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom; %cuidado, ahora está sin capacidad
            end
        end


        for i=1:n
            eta*sum(a(:,i)) <= s(i);
            eta*sum(a(i,:)) <= s(i);
        end
        sum(sum(a)) <= a_max;
        % for o=1:n
        %     for d=1:n
        %         sum(fij(o,:,o,d)) == z(o,d);
        %     end
        % end
        
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

    %tengo que poner las restricciones de la logit, sacar una funcion para
    %obtener las pendientes en funcion de omega y text. probablemente el dm
    %cambie para el problema linealizado, estudiar si lo obtengo o resuelvo
    %para alfa = 1.
    cvx_end


    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;

    %disp(alfa_rod(:,1,2)); disp(delta_rod(:,1,2));
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    opt_gap = 100.*(1e6.*cvx_optbnd-obj_val)./obj_val;
    %disp(f); disp(a);
    filename = sprintf('./results_paper/sol_prueba_MIP_4h_6node_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget','opt_gap');

end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj,z,delta_rod,alfa_rod] = compute_sim_MIP_24h(lam,nom_bud,alfa,dm_pax,dm_op,beta)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_6node_network();

    %esto luego hay que quitarlo, lo de alfa = 1
    %alfa = 1;

    M = 1e4;
    nreg = 10;
    %nreg = 7;

    eps = 1e-3;
    %vals_regs = [0.75,0.5,0.25];
    vals_regs = [0.95,0.9,0.87,0.85,0.83,0.8,0.77,0.75,0.73,0.7,0.67,0.65,0.6,0.55,0.5,0.4,0.3,0.2,0.1];
    
    %vals_regs = linspace(0.95,0.05,6);
    vals_regs = linspace(0.97,0.02,nreg-1);
    [lin_coef,dmax,dmin] = get_linearization(n,nreg,alt_time,alt_price,-1,-1,vals_regs);
    tic;
    
    cvx_begin 
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
        bg = bg + sum(station_capacity_slope'.*s_prim);
        bg = bg + sum(sum(link_capacity_slope.*a_prim));
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                tau.*squeeze(sum(sum(squeeze(permute(zij(i,j,:,:),[3 4 1 2]).*demand)))) <= a(i,j).*a_nom; %cuidado, ahora está sin capacidad
            end
        end


        for i=1:n
            eta*sum(a(:,i)) <= s(i);
            eta*sum(a(i,:)) <= s(i);
        end
        sum(sum(a)) <= a_max;
        % for o=1:n
        %     for d=1:n
        %         sum(fij(o,:,o,d)) == z(o,d);
        %     end
        % end
        
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

    %tengo que poner las restricciones de la logit, sacar una funcion para
    %obtener las pendientes en funcion de omega y text. probablemente el dm
    %cambie para el problema linealizado, estudiar si lo obtengo o resuelvo
    %para alfa = 1.
    cvx_end


    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;

    %disp(alfa_rod(:,1,2)); disp(delta_rod(:,1,2));
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    opt_gap = 100.*(1e6.*cvx_optbnd-obj_val)./obj_val;
    %disp(f); disp(a);
    filename = sprintf('./results_paper/sol_prueba_MIP_24h_6node_beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget','opt_gap');

end

function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_6node_network()

    n = 6;
    
    %candidates to construct a link for each neighbor
    candidates = {[2,3],[1,3,4],[1,2,4,5],[2,3,5,6],[3,4,6],[4,5]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.6,0.8,2,1.6,2.5; 
                2,0,0.9,1.2,1.5,2.5; 
                1.5,1.4,0,1.3,0.9,2; 
                1.9,2,1.9,0,1.8,2; 
                3,1.5,2,2,0,1.5; 
                2.1,2.7,2.2,1,1.5,0];
    
    %fixed cost for constructing links
    link_cost = (1e6/(25*365.25)).*[0,1.7,2.7,0,0,0; 
                 1.7,0,2.1,3,0,0;
                 2.7,2.1,0,2.6,1.7,0; 
                 0,3,2.6,0,2.8,2.4; 
                 0,0,1.7,2.8,0,1.9; 
                 0,0,0,2.4,1.9,0];
    link_cost (link_cost ==0) = 1e4;

    
    %fixed cost for constructing stations
    station_cost = (1e6/(25*365.25)).*[2, 3, 2.2, 3, 2.5, 1.3];
    
    link_capacity_slope = 0.04.*link_cost; 
    station_capacity_slope = 0.04.*station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,9,26,19,13,12;
              11,0,14,26,7,18;
              30,19,0,30,24,8;
              21,9,11,0,22,16;
              14,14,8,9,0,20;
              26,1,22,24,13,0];
    
    distance = 10000 * ones(n, n); % Distances between arcs
    
    for i = 1:n
        distance(i, i) = 0;
    end
    
    distance(1, 2) = 0.75;
    distance(1, 3) = 0.7;

    
    distance(2, 3) = 0.6;
    distance(2, 4) = 1.1;
    
    distance(3, 4) = 1.1;
    distance(3, 5) = 0.5;

    
    distance(4,5) = 0.8;
    distance(4,6) = 0.7;

    
    distance(5,6) = 0.5;

    
    for i = 1:n
        for j = i+1:n
            distance(j, i) = distance(i, j); % Distances are symmetric
        end
    end
    
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