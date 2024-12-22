close all;
clear all;
clc;
%poner logit coef a 1/10
[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
    a_max,candidates,pond_coefs] = parameters_sevilla_network();

%% Resolvemos los problemas por separado para calcular el punto de utopía

%Problema del operador:

%op_obj = operation_costs
%obj_val = 0


%Problema de los pasajeros:
%pax_obj = (distances + prices) + entropies + congestion
%obj_val. Se calcula para un presupuesto dado. En este caso daremos un
%valor a beta, calculamos el presupuest obtenido y lo mantenemos todo el
%tiempo.
betas = [0,0.01,0.02,0.04,0.06,0.08,0.1,1,10,100,1e3,1e4];
%lam 20: entre 0.08 y 0.1
%lam 10: entre 0.1 y 1. 0.12,0.14,0.16,0.18
%lam 7: entre 0.1 y 1. 0.15,  0.22,0.24,0.26,0.28
%lam 5: entre 0.1 y 1. 0.22,0.24,0.26,0.28
%lam 3: entre 0.1 y 1. 0.35, 0.41, 0.43, 0.45
%lam 1: entre 0.1 y 1. quiza alguno > 1. 0.75 y 0.95
betas = 0.2:0.1:0.9;
betas = [0.22,0.24,0.26,0.28];
betas = [0.35,0.41,0.43,0.45,0.75,0.95];
betas = [0.12,0.15,0.23,0.35,0.41,0.43,0.45];
betas = [0.17];
betas = [0.01,0.1,0.12,0.15,0.2,0.22,0.23,0.24,0.26,0.28,0.3,0.35,0.4,0.41,0.43,0.45,0.5,0.6,0.7,0.8,0.9,1,10];

betas = [0.35,0.4,0.41,0.43,0.45,0.5,0.6,0.7,0.8,0.9,1];

betas = [0.01,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7];
betas = [0.01,0.1,0.15,0.2,0.25,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7];
betas = [0.3];
betas = [0.57,0.425,0.32,0.26,0.27,0.22,0.23];
betas = [0.01,0.1,0.15,0.2,0.22,0.23,0.25,0.26,0.27,0.3,0.32,0.35,0.4,0.425,0.45,0.5,0.55,0.57,0.6,0.65,0.7];
betas = [0.12,0.13,0.17,0.21,0.28];
betas = [0.105,0.11,0.14,0.16];
betas = [0.113,0.115];
lams = [5,7,10,12,15];
betas = [0.01,0.1,0.12,0.15,0.2,0.22,0.23,0.24,0.26,0.28,0.3,0.35,0.4,0.41,0.43,0.45,0.5,0.6,0.7,0.8,0.9,1,10];
betas = [0.125,0.13,0.16,0.17,0.32,0.33];
betas = [0.14,0.18,0.19];
betas = [0.132];


lams = [7,10,12,15];
lams = [15];
betas = [0.11,0.15,0.2,0.25];
betas = [0.12,0.125,0.13];
betas = [0.1,0.09,0.11];
betas = 0.06:0.01:0.08;
betas = 0.09:0.01:0.11;
%betas = 0.082:0.02:0.088;
betas = 0.084:0.002:0.088;
betas = 0.092:0.002:0.098;
betas = 0.0905:0.0005:0.0915;

%these are the ones for the Sevilla simulations in the paper for lam = 15
lams = [15];
betas = 0.05:0.01:0.1;
betas = [0.03,0.04,0.073,0.076,0.081];
betas = 0.082:0.001:0.089;
% low income passengers priorization
betas = 0.02:0.01:0.12;
betas = 0.13:0.01:0.25;
betas = [0.161:0.001:0.169,0.171:0.001:0.179];


%these are the ones for the Sevilla simulations in the paper for lam = 10
 lams = [10];
 betas = 0.05:0.01:0.15;
 betas = 0.12:0.002:0.14;
 betas = [0.125,0.127,0.1305,0.131,0.1315];
 betas = [0.1272,0.1275,0.1277];
 % low income passengers priorization
 betas = 0.02:0.01:0.12;
 betas = 0.13:0.01:0.25;
 betas = [0.232,0.234,0.236,0.238,0.252:0.002:0.27];
 betas = [0.2602,0.2604,0.2606,0.2608,0.261,0.2612,0.2614,0.2616,0.2618,0.2625,0.263,0.2635,0.265]
% 
% %these are the ones for the Sevilla simulations in the paper for lam = 5
 lams = [5];
  % low income passengers priorization
 betas = 0.1:0.02:0.3;
 betas = 0.32:0.02:0.5;
 betas = [0.385,0.39,0.395,0.445,0.45,0.455,0.462:0.002:0.478];
 betas = [0.4475,0.456,0.457,0.458,0.459,0.4605,0.461,0.4615,0.4645,0.465,0.4655,0.467,0.471];
% betas = 0.1:0.01:0.2;
% betas = 0.2:0.01:0.4;
% betas = [0.215,0.221,0.222,0.223,0.224,0.225,0.226,0.227,0.228,0.229,...
%     0.231,0.232,0.233,0.234,0.235,0.236,0.237,0.238,0.239];

%betas = [0.3];
%lams = [1,3,5,7,10,20];


cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision high
cvx_save_prefs 

%num_workers = 8;
% Inicia el parpool (piscina de trabajadores) para parfor
% parpool('local', num_workers); % num_workers es el número de trabajadores a utilizar
% parpool(num_workers);
lam = lams(1);

%parfor ll=1:length(lams)
 %   lam = lams(ll);
    parfor bb=1:length(betas)
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
delete(gcp);

%% See the results
clear all;
betas = 0.2:0.1:0.9;
betas = [0.22,0.24,0.26,0.28];
betas = [0.35,0.41,0.43,0.45,0.75,0.95];
betas = [0,0.01,0.02,0.04,0.06,0.08,0.1,1,10,100,1e3,1e4];
%betas = [0.01,0.1,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.41,0.43,0.45,0.5,0.6,0.7,0.75,0.8,0.9,0.95,1];
%with poverty
betas = [0.01,0.1,0.15,0.2,0.22,0.23,0.25,0.26,0.27,0.3,0.32,0.35,0.4,0.425,0.45,0.5,0.55,0.57,0.6,0.65,0.7];
betas = [0.15];
%betas = [0.1];
%poverty:
%lam 5: 0.57
%lam 7: 0.425
%lam 10: 0.32
%lam 12: 0.26,0.27
%lam 15: 0.22,0.23

%original:
betas = [0.01,0.1,0.105,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.2,0.21,0.22,0.23,0.25,0.26,0.27,0.28,0.3,0.32,0.35,0.4,0.425,0.45,0.5,0.55,0.57,0.6,0.65,0.7];
 %betas = [0.13];
%lam 5: 0.27-0.3
%lam 7: 0.2-0.22
%lam 10: 0.15- 0.16
%lam 12: 0.14
%lam 15: 0.11, 0.105

% betas = [0.01,0.1,0.12,0.125,0.13,0.132,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.22,0.23,0.24,0.26,0.28,0.3,0.32,0.33,0.35,0.4,0.41,0.43,0.45,0.5,0.6,0.7,0.8,0.9,1];
%
%betas = [0.16];
%crecimiento:
%lam 5: 
%lam 7: 
%lam 10: 
%lam 12: 
%lam 15: 0.132
%betas = [0.132];


lam = 10;
close all;

betas = [0.15];

lams = [7,10,12,15];
lams = [7];
betas = [0.11,0.15,0.2,0.25];
betas = [0.11,0.12,0.125,0.13];
betas = [0.1,0.09,0.11];
betas = 0.06:0.01:0.08;
%
betas = 0.05:0.01:0.1;
betas = [0.03,0.04,0.073,0.076,0.081:0.001:0.089];
betas = [0.073];
%betas = [0.0915];
%betas = [0.11];
lams = [15];
lams = [10];
betas = 0.05:0.01:0.15;
betas = 0.12:0.002:0.14;
betas = [0.125,0.127,0.1305,0.131,0.1315];
betas = [0.1272,0.1275,0.1277];
betas = 0.1:0.01:0.2;
betas = 0.2:0.01:0.4;
betas = [0.215,0.221,0.222,0.223,0.224,0.225,0.226,0.227,0.228,0.229,...
    0.231,0.232,0.233,0.234,0.235,0.236,0.237,0.238,0.239];
%0:22:0.001:0.24;
betas = 0.02:0.01:0.12;
betas = 0.13:0.01:0.25;
betas = [0.161:0.001:0.169,0.171:0.001:0.179];
 betas = 0.02:0.01:0.12;
  betas = 0.13:0.01:0.25;
   betas = [0.232,0.234,0.236,0.238,0.252:0.002:0.27];
    betas = [0.2602,0.2604,0.2606,0.2608,0.261,0.2612,0.2614,0.2616,0.2618,0.2625,0.263,0.2635,0.265]
lam = 10;

lam = 5;
 betas = 0.1:0.02:0.3;
  betas = 0.32:0.02:0.5;
   betas = [0.385,0.39,0.395,0.445,0.45,0.455,0.462:0.002:0.478];
 betas = [0.4475,0.456,0.457,0.458,0.459,0.4605,0.461,0.4615,0.4645,0.465,0.4655,0.467,0.471];

%baseline
lam = 5; betas = 0.21;
lam = 10; betas = 0.12;
lam = 15; betas = 0.08;

lam=5; betas= 0.224;
lam=10; betas = 0.128;
lam = 10; betas = 0.14;
% lam=15; betas = 0.085;
% %pobreza
% lam = 5; betas = 0.4;
% lam = 10; betas = 0.232;
% lam=15; betas = 0.163;
% 
% lam = 5; betas = 0.464;
% lam = 10; betas = 0.262;
% lam = 15; betas = 0.172;

clc;
for bb=1:length(betas)
    beta_or = betas(bb);
    [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates] = parameters_sevilla_network();

    filename = sprintf('./results_paper/original_centro_caro_river_beta=%d_lam=%d.mat',beta_or,lam);
   % filename = sprintf('./sevilla_lines_trial/original_centro_caro_poverty_river_beta=%d_lam=%d.mat',beta_or,lam);

    load(filename);
    att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
    bud(bb) = budget;
    nroutes = 0;
    dis_rut = 0;
    for o=1:n
        for d=1:n
            dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
            if f(o,d) > 0.01
                dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                nroutes = nroutes + 1;
            end
            uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
        end
    end
    long_mean(bb) = dis_rut/nroutes;
    d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
    u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
    disp(['gamma = ',num2str(1/beta_or)]);
    disp(['Dem att. = ',num2str(100*sum(sum(f.*demand))/sum(sum(demand))),' %']);
    att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
   % disp(['att_dem = ',num2str(sum(sum(f.*demand))),' pax']);
    disp(['Pres. = ',num2str(budget)]);
    disp(['Arcos = ',num2str(sum(sum(a_prim > 0.1)))]);
    disp(['Nodos = ',num2str(sum(s_prim > 0.1))]);
    disp(['Cap. = ',num2str(sum(sum(a_prim)))]);
    disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
    disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
    disp(['distancia por ruta = ',num2str(long_mean(bb))]);
    disp(['tiempo computacional = ',num2str(comp_time)]);
    disp('\n');
    
end

figure(1);
subplot(211);
plot(bud,att_d,'-o','LineWidth',1.5); grid on; yl = ylabel('demanda atraída [%]'); xl = xlabel('presupuesto [€]');
xl.FontSize = 15;
yl.FontSize = 15;
title("Demanda atraída VS Presupuesto",'FontSize',15);
subplot(212);
plot(1./betas,att_d,'-o','LineWidth',1.5); grid on; xl = xlabel('\gamma'); yl = ylabel('demanda atraída [%]');
xl.FontSize = 15;
yl.FontSize = 15;
title("Demanda atraída VS \gamma",'FontSize',15);
tit = num2str(['\lambda = ',num2str(lam)]);
sgtitle(tit);
figure(2);
subplot(211);
plot(1./betas,d_med,'-o','LineWidth',1.5); 
hold on
plot(1./betas,u_med,'-o','LineWidth',1.5); grid on; xl = xlabel('\gamma'); yl = ylabel('Tiempo medio por pasajero');
legend('Red nueva','Red competidora');
xl.FontSize = 15;
yl.FontSize = 15;
subplot(212);
plot(1./betas,long_mean,'-o','LineWidth',1.5); grid on; xl = xlabel('\gamma'); yl = ylabel('Tiempo medio por ruta');
tit = num2str(['\lambda = ',num2str(lam)]);
xl.FontSize = 15;
yl.FontSize = 15;
sgtitle(tit);



long = [-5.99652514, -5.98410857, -5.98975246, -5.9400862 , -6.02587338,-5.89719261, -6.00104025, -6.00216903, -5.99878269, -6.00555537,...
    -5.92879841, -6.01571438, -6.03377484, -5.97959346, -5.96040422, -5.97620712, -6.02022949, -5.98410857, -6.00310964, -5.97191776,...
    -5.98884944, -5.91807502, -5.96068641, -6.0303885 ];
lat = [37.3525353 , 37.36021721, 37.3790728 , 37.35323365, 37.39234155, 37.42446589, 37.41119714, 37.38116786, 37.38326293, 37.38326293, ...
    37.40700701, 37.41399056, 37.39373826, 37.39443661, 37.38745306, 37.38815141, 37.3930399 , 37.37767609, ...
    37.40107099, 37.41790135,  37.38926878, 37.40962584, 37.36615323, 37.3790728 ];

lat(9) = 37.3820073;
long(9) = -5.9941967;
 long(21) = -5.9959814;
 long(17) = long(17) + 0.004;
 lat(7) = lat(7) - 0.002;

figure(3);        
a_b = zeros(n);
a_b(a>1) = 1;
a_b = a_b + a_b';
a_b (a_b < 1.9) = 0;
g = graph(a_b);
colormap parula
colors = colormap;
est_size = 3*ones(1,n);
for i=1:n
    if s_prim(i) > 0
        est_size(i) = 4.*(s_prim(i).^0.2);
    end
end
colores_edg = colors( round(1e-2 + g.Edges.Weight .* length(colors)/max(g.Edges.Weight)) ,:);
h = plot(g,'XData',long-mean(long),'YData',lat-mean(lat),'MarkerSize',est_size,'LineWidth',0.7.*(g.Edges.Weight).^0.7,'EdgeColor','blue','EdgeAlpha',1);

pond_coefs = [1.54535294, 1.54535294, 1.11218247, 1.72777031, 2.75490793, ...
        2.14294615, 1.41994992, 1, 1.11218247, 1.41994992,...
           2.14294615, 2.35701223, 2.75490793, 2.21242079, 3.        ,...
           1.59083705, 1.41994992, 1.01231062, 1.41994992, 2.35701223,...
           1.11218247, 2.14294615, 1.72777031, 1.45204153]';

crec_coefs = [1.6, 1.6, 1.1469802107427398, 0.9, 1.2081587290019313, 0.9661783579052806, 1.0802156586966714, ...
        0.9, 1.1469802107427398, 1.0802156586966714, 0.9661783579052806, 0.9989307690507944,...
        1.2081587290019313, 0.9, 1.0730507219686063, 0.9, 1.0802156586966714, 1.3414585127763425, ...
        1.0802156586966714, 0.9989307690507944, 1.1469802107427398, 0.9661783579052806, 0.9, 1.1224800389355853];
% g = graph(zeros(n));
% h = plot(g,'XData',long-mean(long),'YData',lat-mean(lat),'MarkerSize',7);


%colorbar;
% caxis([min(g.Edges.Weight),max(g.Edges.Weight)]);
hold on
I = imread('./sevilla.png'); 
%I = rgb2gray(I);
h = image(xlim+0.008,-0.8*ylim-0.001,I); 
uistack(h,'bottom');
set(h,'AlphaData',0.8);
xticks([]); yticks([]); 
set(gcf, 'Position', [100, 100, 935, 587]);
%set(gcf, 'Position', [100, 80, 935, 900]);
axis image
axis off
tit = ['Diseño con predicción demográfica. \gamma = ' num2str(round(1/beta_or,1)), ' ,  \lambda = ', num2str(lam), '. Demanda atendida = ',num2str(att_dem),' %'];
%tit = ['Coeficientes asociados a la predicción de demanda'];
%title(tit,'FontSize',15);
figurename = sprintf('./figures/topologia_original_lam=%d.png',round(lam));
%saveas(gcf, figurename);
%tit = ['\lambda = ',num2str(lam),', \beta = ',num2str(beta_or),', ',num2str(sum(s_h > 0)),' hubs',...
%    ', att dem = ',num2str(att_dem*100),'%',' , budg = ',num2str(budget)];
%title(tit);


%%
clear all; close all; clc;

%baseline
lam = 5;

betas = [0.2, 0.21, 0.215, 0.22, 0.221, 0.222, 0.223, 0.224, 0.225, 0.226,...
    0.227, 0.228, 0.229, 0.23, 0.231, 0.232, 0.233, 0.234, 0.235, 0.236,...
    0.237, 0.238, 0.239, 0.24, 0.25];

lam = 10;
betas = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.122, 0.124,...
    0.125, 0.126, 0.127, 0.1272, 0.1275, 0.1277, 0.128, 0.13, 0.1305,...
    0.131, 0.1315, 0.132, 0.14];

lam = 15;
betas = [0.03, 0.04, 0.05, 0.06, 0.07, 0.073, 0.076, 0.08, 0.081,...
    0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.1];

%poverty
lam = 5;
betas = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32,...
    0.34, 0.36, 0.38, 0.385, 0.39, 0.395, 0.4, 0.42, 0.44, 0.445, 0.4475,...
    0.45, 0.455, 0.456, 0.457, 0.458, 0.459, 0.46, 0.4605, 0.461, 0.4615,...
    0.462, 0.464, 0.4645, 0.465, 0.4655, 0.466, 0.467, 0.468, 0.47, 0.471,...
    0.472, 0.474, 0.476, 0.478, 0.48, 0.5];
% 
lam = 10;
betas = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,...
    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.232,...
    0.234, 0.236, 0.238, 0.24, 0.252, 0.254, 0.256, 0.258, 0.26, 0.2602,...
    0.2604, 0.2606, 0.2608, 0.261, 0.2612, 0.2614, 0.2616, 0.2618, 0.262,...
    0.2625, 0.263, 0.2635, 0.264, 0.265, 0.266, 0.268, 0.27];
% 
lam=15;
betas = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,...
    0.12, 0.13, 0.14, 0.15, 0.16, 0.161, 0.162, 0.163, 0.164, 0.165, 0.166,...
    0.167, 0.168, 0.169, 0.17, 0.171, 0.172, 0.173, 0.174, 0.175, 0.176,...
    0.177, 0.178, 0.179, 0.18, 0.19];



clc;
for bb=1:length(betas)
    beta_or = betas(bb);
    [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates] = parameters_sevilla_network();
   % filename = sprintf('./sevilla_lines_trial/original_centro_caro_river_beta=%d_lam=%d.mat',beta_or,lam);
    filename = sprintf('./results_paper/original_centro_caro_poverty_river_beta=%d_lam=%d.mat',beta_or,lam);
    load(filename);
    att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
    bud(bb) = budget;
    nroutes = 0;
    dis_rut = 0;
    for o=1:n
        for d=1:n
            dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
            if f(o,d) > 0.01
                dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                nroutes = nroutes + 1;
            end
            uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
        end
    end
    long_mean(bb) = dis_rut/nroutes;
    d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
    u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
    % disp(['gamma = ',num2str(1/beta_or)]);
    % disp(['Dem att. = ',num2str(100*sum(sum(f.*demand))/sum(sum(demand))),' %']);
    att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
    n_links = sum(sum(a_prim));
    n_nodes = sum(s_prim > 0.1);
    served_demand = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
    gamma = 1/beta_or;
    %served_demand = 0;
   % disp(['att_dem = ',num2str(sum(sum(f.*demand))),' pax']);
    % disp(['Pres. = ',num2str(budget)]);
    % disp(['Arcos = ',num2str(sum(sum(a_prim > 0.1)))]);
    % disp(['Nodos = ',num2str(sum(s_prim > 0.1))]);
    % disp(['Cap. = ',num2str(sum(sum(a_prim)))]);
    % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
    % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
    % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
    % disp(['tiempo computacional = ',num2str(comp_time)]);
    % disp('\n');
    disp([sprintf('%.2f',gamma),'&', ...
        num2str(served_demand), ...
        '&',sprintf('%.3e',budget),'&',num2str(sum(sum(a_prim > 0.1))), ...
        '&',num2str(n_nodes), ...
        '&',num2str(n_links), ...
        '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
        sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);


    
end




    


%%

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
    a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();

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
    n = 24;
    logit_coef = 0.2;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;
    op_obj = op_obj + 1e-3*(sum(sum(op_link_cost.*a_prim))); %operational costs
    for i=1:n
        if s_prim(i) > eps
            pax_obj = pax_obj - 1e-3*sum(log(congestion_coef_stations(i)*delta_s(i) + eps));
        end

        for j=1:n
            if a_prim(i,j) > eps
                pax_obj = pax_obj - 1e-3*sum(sum(log(congestion_coef_links(i,j)*delta_a(i,j) + eps)));
            end
        end
    end
    for o=1:n
        for d=1:n
            for i=1:n
                for j=1:n
                    pax_obj = pax_obj + 1e-3*(demand(o,d).*logit_coef*(travel_time(i,j)+prices(i,j)).*fij(i,j,o,d));
                end
            end
        end
    end
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*logit_coef*(alt_time+alt_price).*fext)));
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(-entr(f) - f))));
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(-entr(fext) - fext))));
    obj_val = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim(lam,beta_or,alfa,dm_pax,dm_op)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
    a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();

    %pond_coefs_mat = ones(n);

    niters = 10;           
    eps = 1e-3;
    a_prev = 1e4.*ones(n);
    s_prev= 1e4.*ones(n,1);
    logit_coef = 0.2;
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
            bud_obj = bud_obj + 1e-3*sum(station_capacity_slope.*s_prim);
            bud_obj = bud_obj + 1e-3*(sum(sum(link_capacity_slope.*a_prim)));  %linear construction costs
            op_obj = op_obj + 1e-3*(sum(sum(op_link_cost.*a_prim))); %operation costs
            if iter < niters
                for i=1:n
                    pax_obj = pax_obj + 1e-3*(sum(inv_pos(congestion_coef_stations(i).*delta_s(i) + eps)));
                    for j=1:n
                        pax_obj = pax_obj + 1e-3*(sum(sum(inv_pos(congestion_coef_links(i,j).*delta_a(i,j) + eps))));
                    end
                end
                bud_obj = bud_obj + 1e-3*lam*sum(sum((link_cost.*a_prim.*(1./(a_prev+eps))))) + 1e-3*lam*sum((station_cost.*s_prim.*(1./(s_prev+eps)))); %fixed construction costs
                
            end
    
            for o=1:n
                for d=1:n
                    pax_obj = pax_obj + 1e-3*(pond_coefs_mat(o,d)*demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d).*logit_coef))); 
                end
            end
            pax_obj = pax_obj + 1e-3*(sum(sum(pond_coefs_mat.*demand.*(alt_time+alt_price).*fext.*logit_coef)));
            pax_obj = pax_obj + 1e-3*(sum(sum(pond_coefs_mat.*demand.*(-entr(f) - f))));
            pax_obj = pax_obj + 1e-3*(sum(sum(pond_coefs_mat.*demand.*(-entr(fext) - fext))));
    
    
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
            obj = beta_or*bud_obj + (alfa/(dm_pax))*pax_obj + ((alfa)/(dm_op))*op_obj;
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

            for i=1:n
                for j=1:n
                    a_prim(i,j) == a_prim(j,i);
                end
            end
    
            for i=1:n
                for j=1:n
                    squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= tau.*a(i,j).*a_nom; %multiplicar por demanda
                end
            end
            for i=1:n
              %  sigma*sum(a_prim(:,i)) <= s(i); %cambiar esta por solo los que salen o entran en la estacion
                sigma*sum(squeeze(permute(sum(fij(i,:,i,:),2),[3 4 1 2])).*demand(i,:))./(tau*a_nom)  + sigma*sum(delta_a(i,:)) <= s(i);
                sigma*sum(squeeze(permute(sum(fij(:,i,:,i),1),[3 4 1 2])).*demand(:,i))./(tau*a_nom) + sigma*sum(delta_a(:,i)) <= s(i);

            end
            sum(sum(a_prim)) <= a_max;
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
                        a_prim(i,j) == 0; %ver esto
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

        disp(['iter = ',num2str(iter),', beta = ',num2str(beta_or),', lam = ',num2str(lam),', nlinks =',num2str(sum(sum(a_prim > 0.1))),', nstations = ', ...
            num2str(sum(s_prim > 0.1)),', att_dem = ',num2str(100*sum(sum(f.*demand))/sum(sum(demand))),' %']);


        if beta_or == 0.08


            long = [-5.99652514, -5.98410857, -5.98975246, -5.9400862 , -6.02587338,-5.89719261, -6.00104025, -6.00216903, -5.99878269, -6.00555537,...
                -5.92879841, -6.01571438, -6.03377484, -5.97959346, -5.96040422, -5.97620712, -6.02022949, -5.98410857, -6.00310964, -5.97191776,...
                -5.98884944, -5.91807502, -5.96068641, -6.0303885 ];
            lat = [37.3525353 , 37.36021721, 37.3790728 , 37.35323365, 37.39234155, 37.42446589, 37.41119714, 37.38116786, 37.38326293, 37.38326293, ...
                37.40700701, 37.41399056, 37.39373826, 37.39443661, 37.38745306, 37.38815141, 37.3930399 , 37.37767609, ...
                37.40107099, 37.41790135,  37.38926878, 37.40962584, 37.36615323, 37.3790728 ];
            
            lat(9) = 37.3820073;
            long(9) = -5.9941967;
             long(21) = -5.9959814;
             long(17) = long(17) + 0.004;
             lat(7) = lat(7) - 0.002;
            
            figure(3);        
            a_b = zeros(n);
            a_b(a>1) = 1;
            a_b = a_b + a_b';
            a_b (a_b < 1.9) = 0;
            g = graph(a_b);
            colormap parula
            colors = colormap;
            est_size = 3*ones(1,n);
            for i=1:n
                if s_prim(i) > 0
                    est_size(i) = 4.*(s_prim(i).^0.2);
                end
            end
            colores_edg = colors( round(1e-2 + g.Edges.Weight .* length(colors)/max(g.Edges.Weight)) ,:);
            h = plot(g,'XData',long-mean(long),'YData',lat-mean(lat),'MarkerSize',est_size,'LineWidth',0.7.*(g.Edges.Weight).^0.7,'EdgeColor',colores_edg,'EdgeAlpha',1);
            
            
            %colorbar;
            % caxis([min(g.Edges.Weight),max(g.Edges.Weight)]);
            hold on
            I = imread('./sevilla.png'); 
            %I = rgb2gray(I);
            h = image(xlim+0.008,-0.8*ylim-0.001,I); 
            uistack(h,'bottom');
            set(h,'AlphaData',0.8);
            xticks([]); yticks([]); 
            set(gcf, 'Position', [100, 100, 935, 587]);
            %set(gcf, 'Position', [100, 80, 935, 900]);
            axis image
            axis off

        end

    end
    comp_time = toc;
    
    
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./results_paper/original_centro_caro_beta=%d_lam=%d.mat',beta_or,lam);
    filename = sprintf('./results_paper/original_centro_caro_river_logit_soft_beta=%d_lam=%d.mat',beta_or,lam);
    filename = sprintf('./results_paper/original_centro_caro_river_beta=%d_lam=%d.mat',beta_or,lam);
    filename = sprintf('./results_paper/original_centro_caro_poverty_river_beta=%d_lam=%d.mat',beta_or,lam);
    %filename = sprintf('./sevilla_lines_trial/original_beta=%d_lam=%d.mat',beta_or,lam);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

end


function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
    a_max,candidates,pond_coefs_mat] = parameters_sevilla_network()

    n = 24;
    
    %candidates to construct a link for each neighbor
    candidates = {
    [4,23,2,18,3,9,8,10,24,5,17];
    [1,4,23,18,3,9,8,5,24];
    [1,2,18,16,14,21,19,9,8,24];
    [1,2,23,16,15,20,11,22,6];
    [1,2,10,9,21,14,17,12,13,24];
    [4,22,16,14,19,7,20];
    [12,13,17,19,21,14,11,22,6,20];
    [1,2,3,9,10,24];
    [1,2,3,16,21,20,19,17,5,10,8,12];
    [1,8,9,20,19,12,17,5,24,21];
    [4,22,20,7,19,14,16,15,23];
    [13,5,17,10,9,19,14,7,20];
    [24,5,17,19,7,12];
    [21,3,18,16,11,22,6,20,7,12,19,17,5,15];
    [18,23,4,22,11,20,14,16];
    [18,23,4,15,11,6,14,21,9,3,22];
    [5,24,1,10,9,21,14,19,7,12,13];
    [2,23,15,16,14,21,19,3,1];
    [13,17,10,9,3,18,21,14,11,6,20,7,12];
    [6,22,11,4,15,14,21,9,10,19,7,12];
    [3,18,23,16,14,20,7,19,17,5,10,9];
    [4,6,20,7,14,16,11,15,23];
    [1,4,22,11,15,16,21,18,24,2];
    [1,2,23,3,8,10,17,5,13];
    };
    

    population_file = readtable('./population.xlsx');
    population = table2array(population_file(1:24,2));
    coordinates = readtable('./coordenadas_Sevilla.xlsx');
    coor_x = table2array(coordinates(1:24,3));
    coor_y = table2array(coordinates(1:24,7));
    rng(1,"twister"); %seed
    distance = 1e6.*ones(n);
    for i=1:n
        distance(i,i) = 0;
        cand = candidates(i);
        cand = cand{1};
        cand = cand(cand > i);
        for j=i+1:n
            if sum(j == cand) > 0
                distance(i,j) = haversine(coor_y(i), coor_x(i), coor_y(j), coor_x(j));
                distance(j,i) = distance(i,j);
            end
            non_stop = rand < 0.4;

            alt_cost(i,j) = haversine(coor_y(i), coor_x(i), coor_y(j), coor_x(j));
            alt_cost(i,j) = alt_cost(i,j) + 0.2*non_stop*alt_cost(i,j);
            alt_cost(j,i) = alt_cost(i,j);
        end
    end
    demand = [0, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           272, 0, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           327, 327, 0, 327, 327, 664, 327, 327, 327, 664, 664, 327, 664, 327, 664, 664, 664, 327, 1125, 1125, 1125, 1125, 1125, 1125;
           185, 185, 185, 0, 185, 376, 185, 185, 185, 376, 376, 185, 376, 185, 376, 376, 376, 185, 637, 637, 637, 637, 637, 637;
           272, 272, 272, 272, 0, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           225, 225, 225, 225, 225, 0, 225, 225, 225, 188, 188, 225, 188, 225, 188, 188, 188, 225, 284, 284, 284, 284, 284, 284;
           283, 283, 283, 283, 283, 575, 0, 283, 283, 575, 575, 283, 575, 283, 575, 575, 575, 283, 975, 975, 975, 975, 975, 975;
           272, 272, 272, 272, 272, 553, 272, 0, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           272, 272, 272, 272, 272, 553, 272, 272, 0, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           511, 511, 511, 511, 511, 428, 511, 511, 511, 0, 428, 511, 428, 511, 428, 428, 428, 511, 645, 645, 645, 645, 645, 645;
           225, 225, 225, 225, 225, 188, 225, 225, 225, 188, 0, 225, 188, 225, 188, 188, 188, 225, 284, 284, 284, 284, 284, 284;
           272, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 0, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           306, 306, 306, 306, 306, 257, 306, 306, 306, 257, 306, 257, 0, 257, 306, 306, 306, 257, 387, 387, 387, 387, 387, 387;
           294, 294, 294, 294, 294, 597, 294, 294, 294, 597, 597, 294, 597, 0, 597, 597, 597, 297, 1012, 1012, 1012, 1012, 1012, 1012;
           409, 409, 409, 409, 409, 342, 409, 409, 409, 342, 342, 409, 342, 409, 0, 342, 342, 409, 516, 516, 516, 516, 516, 516;
           511, 511, 511, 511, 511, 428, 511, 511, 511, 428, 428, 511, 428, 511, 428, 0, 428, 511, 645, 645, 645, 645, 645, 645;
           429, 429, 429, 429, 429, 360, 429, 429, 429, 306, 360, 429, 360, 429, 360, 360, 0, 429, 542, 542, 542, 542, 542, 542;
           272, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 0, 937, 937, 937, 937, 937, 937;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 660, 660;
           879, 879, 879, 879, 879, 952, 879, 879, 879, 952, 952, 879, 952, 879, 952, 952, 952, 879, 860, 0, 860, 860, 860, 860;
           715, 715, 715, 715, 715, 775, 715, 715, 715, 775, 775, 715, 775, 715, 775, 775, 775, 715, 700, 700, 700, 700, 700, 700;
           511, 511, 511, 511, 511, 553, 511, 511, 511, 553, 553, 511, 553, 511, 553, 553, 553, 511, 500, 500, 500, 0, 500, 500;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 0, 660;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 660, 0];

    pond_coefs = [1.54535294, 1.54535294, 1.11218247, 1.72777031, 2.75490793, ...
        2.14294615, 1.41994992, 1, 1.11218247, 1.41994992,...
           2.14294615, 2.35701223, 2.75490793, 2.21242079, 3.        ,...
           1.59083705, 1.41994992, 1.01231062, 1.41994992, 2.35701223,...
           1.11218247, 2.14294615, 1.72777031, 1.45204153]';
    
    pond_coefs_tens(1,:,:) = pond_coefs.*ones(1,n);
    pond_coefs_tens(2,:,:) = permute(pond_coefs_tens(1,:,:),[1 3 2]);
    pond_coefs_mat = squeeze(permute(max(pond_coefs_tens(1,:,:),pond_coefs_tens(2,:,:)),[2,3,1]));

    crec_coefs = [1.6, 1.6, 1.1469802107427398, 0.9, 1.2081587290019313, 0.9661783579052806, 1.0802156586966714, ...
        0.9, 1.1469802107427398, 1.0802156586966714, 0.9661783579052806, 0.9989307690507944,...
        1.2081587290019313, 0.9, 1.0730507219686063, 0.9, 1.0802156586966714, 1.3414585127763425, ...
        1.0802156586966714, 0.9989307690507944, 1.1469802107427398, 0.9661783579052806, 0.9, 1.1224800389355853];

    % for o=1:n
    %     for d=1:n
    %         demand(o,d) = crec_coefs(o)*crec_coefs(d)*demand(o,d);
    %     end
    % end

   % demand = demand.*pond_coefs_mat;
    %demand = max(demand,demand');
    
    %fixed cost for constructing links

    link_cost = 1e6.*distance./(365.25*25);
    group_1 = [3,8,9,10,14,16,18,21];
    group_2 = [1,2,7,15,17,19];
    group_3 = [4,5,6,11,12,13,20,22,23,24];

    % 1-1: 3* (1.5 el 1, 1 el 2, 0.5 el 3)
    % 1-2: 2.5*
    % 1-3: 2*
    % 2-2: 2*
    % 2-3: 1.5*
    % 3-3: 1*

    for i=1:n
        if ismember(i,group_1)
            ci = 1.5;
        end
        if ismember(i,group_2)
            ci = 1;
        end
        if ismember(i,group_3)
            ci = 0.5;
        end

        for j=i+1:n
            if ismember(j,group_1)
                cj = 1.5;
            end
            if ismember(j,group_2)
                cj = 1;
            end
            if ismember(j,group_3)
                cj = 0.5;
            end
            link_cost(i,j) = link_cost(i,j).*(ci+cj);
            link_cost(j,i) = link_cost(i,j);
        end

    end





    nodos_centricos = [3,8,9,10,14,16,17,18,19,21];
    % for i=1:(length(nodos_centricos)-1)
    %     for j=i+1:length(nodos_centricos)
    %         link_cost(nodos_centricos(i),nodos_centricos(j)) = link_cost(nodos_centricos(i),nodos_centricos(j)) + 1e6.*5./(365.25.*25);
    %         link_cost(nodos_centricos(j),nodos_centricos(i)) =link_cost(nodos_centricos(i),nodos_centricos(j));
    %     end
    % end
    
    %fixed cost for constructing stations
   
  

    river_plus = 100;
   % river_plus = 0;


    link_cost(12,[17,10,9,19,14,7,20]) =  link_cost(12,[17,10,9,19,14,7,20]) + river_plus;
    link_cost([17,10,9,19,14,7,20],12) = link_cost([17,10,9,19,14,7,20],12) + river_plus;

    link_cost(5,[1,2,10,9,21,14,17]) = link_cost(5,[1,2,10,9,21,14,17]) + river_plus;
    link_cost([1,2,10,9,21,14,17],5) = link_cost([1,2,10,9,21,14,17],5) + river_plus;

    link_cost(13,[17,19,7]) = link_cost(13,[17,19,7]) + river_plus;
    link_cost([17,19,7],13) = link_cost([17,19,7],13) + river_plus;

    link_cost(24,[1,2,23,3,8,10,17]) = link_cost(24,[1,2,23,3,8,10,17]) + river_plus;
    link_cost([1,2,23,3,8,10,17],24) = link_cost([1,2,23,3,8,10,17],24) + river_plus;

    link_cost(7,[21,14,11,22,6,20]) = link_cost(7,[21,14,11,22,6,20]) + river_plus;
    link_cost([21,14,11,22,6,20],7) = link_cost([21,14,11,22,6,20],7) + river_plus;

    link_cost(19,[9,3,21,18,14,11]) = link_cost(19,[9,3,21,18,14,11]) + river_plus;
    link_cost([9,3,21,18,14,11],19) = link_cost([9,3,21,18,14,11],19) + river_plus;

    link_cost(17,[9,21,14]) =  link_cost(17,[9,21,14]) + river_plus;
    link_cost([9,21,14],17) =  link_cost([9,21,14],17) + river_plus;

    link_cost(10,[9,20,21]) = link_cost(10,[9,20,21]) + river_plus;
    link_cost([9,20,21],10) = link_cost([9,20,21],10) + river_plus;

    link_cost(8,[2,3,9]) = link_cost(8,[2,3,9]) + river_plus;
    link_cost([2,3,9],8) = link_cost([2,3,9],8) + river_plus;

    link_cost(1,[4,23,2,18,3,9]) = link_cost(1,[4,23,2,18,3,9]) + river_plus;
    link_cost([4,23,2,18,3,9],1) = link_cost([4,23,2,18,3,9],1) + river_plus;

    % link_cost(12,[9,14,20]) = link_cost(12,[9,14,20]) + river_plus;
    % link_cost([9,14,20],12) = link_cost([9,14,20],12) + river_plus;
    % 
    % link_cost(5,[2,9,21,14]) = link_cost(5,[2,9,21,14]) + river_plus;
    % link_cost([2,9,21,14],5) = link_cost([2,9,21,14],5) + river_plus;
    % 
    % link_cost(24,[2,23,3,15]) = link_cost(24,[2,23,3,15]) + river_plus;
    % link_cost([2,23,3,15],24) = link_cost([2,23,3,15],24) + river_plus;

    link_capacity_slope = 0.3.*link_cost; 

    station_cost = 1e3.*population./(365.25*25);
    

    station_capacity_slope = 0.2.*station_cost;
    
    
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
    sigma = 0.25;
    a_max = 1e9;
    eps = 1e-3;
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
