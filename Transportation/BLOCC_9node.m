%code adopted from https://github.com/Liuyuan999/Penalty_Based_Lagrangian_Bilevel

%% Inicializacion red 9 nodos

clear all;
close all;
clc;


n = 9;

distance = 10000 * ones(n, n); % Distances between arcs

for i = 1:n
    distance(i, i) = 0;
end

distance(1, 2) = 0.75;
%distance(1,2) = 0.9;
distance(1, 3) = 0.7;
distance(1, 9) = 0.9;

distance(2, 3) = 0.6;
distance(2, 4) = 1.1;

distance(3, 4) = 1.1;
distance(3, 5) = 0.5;
distance(3, 9) = 0.7;

distance(4,5) = 0.8;
distance(4,6) = 0.7;
distance(4,8) = 0.8;
%distance(4,8) = 1.8;

distance(5,6) = 0.5;
distance(5,7) = 0.7;

distance(6,7) = 0.5;
distance(6,8) = 0.4;

for i = 1:n
    for j = i+1:n
        distance(j, i) = distance(i, j); % Distances are symmetric
    end
end

 d = distance;

 u = [0,1.6,0.8,2,1.6,2.5,3,2.5,0.8; ...
    2,0,0.9,1.2,1.5,2.5,2.7,2.4,1.8; ...
    1.5,1.4,0,1.3,0.9,2,1.6,2.3,0.9; ...
    1.9,2,1.9,0,1.8,2,1.9,1.2,2; ...
    3,1.5,2,2,0,1.5,1.1,1.8,1.7; ...
    2.1,2.7,2.2,1,1.5,0,0.9,0.9,2.9; ...
    2.8,2.3,1.5,1.8,0.9,0.8,0,1.3,2.1;...
    2.8,2.2,2,1.1,1.5,0.8,1.9,0,0.3; ...
    1,1.5,1.1,2.7,1.9,1.8,2.4,3,0];

 c = [0,1.7,2.7,0,0,0,0,0,2.9; ...
     1.7,0,2.1,3,0,0,0,0,0; ...
     2.7,2.1,0,2.6,1.7,0,0,0,2.5; ... 
     0,3,2.6,0,2.8,2.4,0,3.2,0; ...
     0,0,1.7,2.8,0,1.9,3,0,0; ...
     0,0,0,2.4,1.9,0,2.7,2.8,0; ...
     0,0,0,0,3,2.7,0,0,0; ...
     0,0,0,3.2,0,2.8,0,0,0; ...
     2.9,0,2.5,0,0,0,0,0,0];

c(c==0) = 1e1;

m = [0,2.5,2.5,5,5,7.5,7.5,7.5,2.5;...
    2.5,0,2.5,2.5,5,5,7.5,5,5;...
    2.5,2.5,0,2.5,2.5,5,5,7.5,2.5;...
    5,2.5,2.5,0,2.5,2.5,5,2.5,7.5;...
    5,5,2.5,2.5,0,2.5,2.5,5,2.5;...
    7.5,5,5,2.5,2.5,0,2.5,2.5,7.5;...
    7.5,7.5,5,5,2.5,2.5,0,5,7.5;...
    7.5,5,5,2.5,5,2.5,5,0,7.5;...
    2.5,5,2.5,5,5,7.5,7.5,7.5,0];

demand = 1e-1.*[0,9,26,19,13,12,13,8,11;
          11,0,14,26,7,18,3,6,12;
          30,19,0,30,24,8,15,12,5;
          21,9,11,0,22,16,25,21,23;
          14,14,8,9,0,20,16,22,21;
          26,1,22,24,13,0,16,14,12;
          8,6,9,23,6,13,0,11,11;
          9,2,14,20,18,16,11,0,4;
          8,7,11,22,27,17,8,12,0];



%% Inicializacion red Sevilla
%   clear all; close all; clc;
% % % 
%  [n,link_cost,station_cost,link_capacity_slope,...
%      station_capacity_slope,demand,prices,...
%      op_link_cost,congestion_coef_stations,...
%      congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
%      a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
% 
%      coordinates = readtable('../red_sevilla/coordenadas_Sevilla.xlsx');
%     coor_x = table2array(coordinates(1:24,3));
%     coor_y = table2array(coordinates(1:24,7));
% 
% 
% 
% 
% demand = demand./2e2;
% 
% 
% adjacency_matrix = zeros(n, n);
% k = 3;
% 
% for i = 1:n
%     % Ordenar las distancias y obtener los índices correspondientes
%     [~, sorted_indices] = sort(travel_time(i, :));
%     % Conservar solo los k-vecinos más cercanos
%     k_nearest_indices = sorted_indices(2:k+1); % El primer vecino es el propio nodo, por eso se omite
%     % Establecer las conexiones en la matriz de adyacencia
%     adjacency_matrix(i, k_nearest_indices) = 1;
%     adjacency_matrix(k_nearest_indices, i) = 1; % La matriz de adyacencia es simétrica
% end
% 
% adjacency_matrix = (adjacency_matrix + adjacency_matrix')./2;
% adjacency_matrix(adjacency_matrix > 0.5) = 1;
% adjacency_matrix(adjacency_matrix <= 0.5) = 0;
% adjacency_matrix(travel_time > 7) = 0;
% 
% travel_time(adjacency_matrix == 0) = 1e2;
% 
% d = 0.25.*travel_time;
% u = 0.25.*alt_time;
% link_cost(adjacency_matrix == 0) = 1e8;
% c = link_cost./2e2;
% 
% for i=1:n
%     for j=1:n
%         if adjacency_matrix(i,j) > 0
%             adjacency_matrix(i,j) = d(i,j);
%         end
%     end
% end
% 
% 
% g = graph(adjacency_matrix);
% plot(g,'XData',coor_x,'YData',coor_y);
% 
% dis_matrix = zeros(n);
% for i=1:n
%     for j=[1:(i-1),(i+1):n]
%         [~, dis_matrix(i,j)] = shortestpath(g, i, j);
%     end
% end
% 
% rng(1);
% b1 = rand(n);
% b1(b1 > 0.5) = 1;
% b1(b1 <= 0.5) = -1;
% b2 = rand(n);
% b2(b2 > 0.5) = 2;
% b2(b2 <= 0.5) = -2;
% 
% m = dis_matrix + b1 + b2;
% m = max(m,0.25);
% 
% 

%% Check computational times for different initializations for inner problem

gamma = 4;
max_iters = 2e5;

rng(1);


a_start = 3 + 12.*rand(n);
tic;
beta = 1.6e-4;
[a_best, fg_best, fijg_best,lamg_best,mug_best,a_vec,f_vec,fij_vec,lamf_vec,muf_vec,obj_vec,obj_vec_f, his_time, best_time] = A3_ods(n,d,u,gamma,m,c,a_start,beta,demand);



figure;
plot(his_time, obj_vec, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('UL utility');
title('9 Node');
grid on;

% 计算并打印最大值
max_obj = max(obj_vec);
max_obj_f = max(obj_vec_f);
fprintf('max obj_vecs: %.6f\n', max_obj);
fprintf('max obj_vecs_f: %.6f\n', max_obj_f);

disp(best_time);

max_iter = 100;
[f,fij,lam,mu, error_f, error_fij] = ll_blocc(n,d,u,a_best,fg_best, fijg_best,lamg_best,mug_best,demand, max_iter);

% %% Results representation
% 
% 
% close all;
% 
% map_x = [2,6,6,11,11,13,14,14,1];
% map_y = [4,6,2,6,2,4,1,7,1];
% 
% 
% figure(1);
% for i=1:length(x)
%     subplot(2,2,i);
%    % imagesc(a_vecs(:,:,i));
%    a = a_vecs(:,:,i);
%    a(a < 0.01) = 0;
%    a = a  + a';
%    g = graph(a);
%    % h = plot(g);
% 
%     h = plot(g,'XData',map_x,'YData',map_y,'LineWidth',0.7.*g.Edges.Weight.^0.7,'NodeFontSize',5,...
%     'EdgeColor','#0072BD','EdgeAlpha',0.8,'interpreter','latex');%, ...
%          %   'MarkerSize',0.7*(s_h+s).^0.5 +1e-2);%,'LineWidth', ...
%             %0.1.*g.Edges.Weight,'NodeColor',colores,'EdgeColor',colores_edg,'EdgeAlpha',0.7,'NodeFontSize',8);
%         xticks([]); yticks([]); 
% 
%    tit = sprintf('$ \\gamma = %d $, $ A_0 =  %.2f$, obj = $ %.2f$, t = $ %.2f$',x(i),1,max(obj_vecs(i,2:end)),comp_time(i));
%    title(tit,'Interpreter','latex','FontSize',5);
% 
%    % title(tit,'FontSize',8);
% end
% 
% figure(2);
% for i=1:length(x)
%     subplot(2,2,i);
%     imagesc(f_vecs(:,:,i));
%     colorbar;
%     hold on;
%     for ii = 1:n
%         for jj = 1:n
%             text(jj, ii, num2str(f_vecs(ii, jj,i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'white');
%             hold on;
%         end
%     end
% end





% 
% %% Lower level for different configurations of the network
% 
% close all;
% xx = 0:0.01:1; x_x = 1-xx;
% y1 = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,2) + d(2,4) + d(4,8)).*xx;
% xx = fg_vec(1); x_x = 1-xx;
% y1_op = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,2) + d(2,4) + d(4,8)).*xx;
% 
% 
% 
% xx = 0:0.01:1; x_x = 1-xx;
% y2 = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,4) + d(4,8)).*xx;
% xx = fg_vec(2); x_x = 1-xx;
% y2_op = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,4) + d(4,8)).*xx;
% 
% 
% xx = 0:0.01:1; x_x = 1-xx;
% y3 = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,5) + d(5,6) + d(6,8)).*xx;
% xx = fg_vec(3); x_x = 1-xx;
% y3_op = xx.*(log(xx) - 1) + x_x.*(log(x_x) - 1) + u.*(x_x) + (d(1,3) + d(3,5) + d(5,6) + d(6,8)).*xx;
% 
% 
% xx = 0:0.01:1; x_x = 1-xx;
% subplot(311); 
% plot(xx,y1); hold on; scatter(fg_vec(1),y1_op,'red','filled'); xlabel('F'); ylabel('PAX(A,F)'); title('A_0 = 0.1')
% subplot(312); plot(xx,y2); hold on; scatter(fg_vec(2),y2_op,'red','filled'); xlabel('F'); ylabel('PAX(A,F)'); title('A_0 = 0.6')
% subplot(313); plot(xx,y3); hold on; scatter(fg_vec(3),y3_op,'red','filled'); xlabel('F');ylabel('PAX(A,F)'); title('A_0 = 1')



%% Functions - all ods



function [f,fij,lam,mu] = A2_ods(n,d,u,a,f_last,fij_last,lam_last,mu_last,demand)
    lam = lam_last;
    f = f_last;
    fij = fij_last;
    mu = mu_last;

    mu_o = zeros(n);
    mu_d = zeros(n);

    mu_i = zeros(n,n,n,n);
    mu_j = mu_i;

    lam_prev = -10*ones(n); %i,j
    f_prev = -10*ones(n); %od
    fij_prev = -10*ones(n,n,n,n); %i,j,o,d
    mu_prev = -10*ones(n,n,n); %i,o,d

    beta_1 = 1e-2;
    beta_2 = 4e-1;
    q = 0;
    max_dif = 1;

    while ( mean(mean(abs(f-f_prev)))/(beta_1) > 2e-2 ) ...
       || ( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1) > 1e-2 ) ...
       || (cons_f_a > 10) ...
       || (max_dif > 2e-1)
               % (mean(min((squeeze(sum(sum(fij,4),3)) < (a+3e-2)))) < 0.5) ... %sin demanda

        lam_prev = lam;
        f_prev = f;
        fij_prev  = fij;
        mu_prev = mu;
   
     %   disp(max_dif)

       

       mu = mu + beta_2.*(squeeze(sum( fij , 2 ))  - ...
          squeeze(permute(sum(fij , 1 ),[2 1 3 4])));

%         for o=1:n
%             for des=2
%                 mu(:,o,des) = mu(:,o,des) + beta*(sum( fij(:,:,o,des),2 ) - sum(fij(:,:,o,des) , 1)' );
%             end
%         end

        for o=1:n
            for des=1:n
                mu(o,o,des) = mu(o,o,des) + beta_2*(-f(o,des));
                mu_o(o,des) = mu(o,o,des);
                mu(des,o,des) = mu(des,o,des) + beta_2*(f(o,des));
                mu_d(o,des) = mu(des,o,des);
            end
        end

     %   mu = max(mu,0);

        for o=1:n
            for des=1:n
                mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
                mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
            end
        end

        for i=1:n
            for j=1:n
                lam(i,j) = lam(i,j) + beta_2.*(   sum(sum(squeeze(permute(fij(i,j,:,:),[3,4,1,2]).*demand))) - a(i,j) );
            end
        end

%         for i=1:n
%             for j=1:n
%                 lam(i,j) = lam(i,j) + beta*(  sum(sum(squeeze(fij(i,j,:,:)))) - a(i,j));
%             end
%         end

        lam = max(0,lam);

        for i=1:n
            lam(i,i) = 100;
        end

        f = f - beta_1*(demand.*log(f./(1-f))  - demand.*u - mu_o + mu_d );

%         for o=1:n
%             for des=1:n
%                 f(o,des) = f(o,des) - beta*(log(f(o,des)./(1-f(o,des)))  - u(o,des) - mu(o,o,des) + mu(des,o,des) ); %puedo quitar algun bucle?
%             end
%         end

        
        f = min(0.99,f);
        f = max(0.01,f);

        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( demand(o,des).*d + lam + mu_i(:,:,o,des) - mu_j(:,:,o,des) );
            end
        end

%         for o=1:n
%             for des=1:n
%                 for i=1:n
%                     for j=1:n
%                         fij(i,j,o,des) = fij(i,j,o,des) - beta*(d(i,j) + lam(i,j) + mu(i,o,des) - mu(j,o,des)  );
%                     end
%                 end
%                % fij(:,:,o,des) = fij(:,:,o,des) - beta*(d + lam + mu_i - mu_j);
%             end
%         end

        fij = min(1,fij);
        fij = max(0,fij);

        for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
            mu(:,i,i) = 0;
        end

        cons_f_a = 0;
        for i=1:n
            for j=1:n
                if squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) >= (a(i,j) + 2e-1)
                    cons_f_a = cons_f_a + 1;
                end
            end
        end



        max_dif = 0;
        
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd))  );
                end
                
            end
        end


        if q == 100
%             disp('fij 1,2 = ');
%             disp(fij(:,:,1,2));
%             disp('f 1,2 = ');
            %disp(f);
            q = 0;
        else
            q = q+1;
        end
        
        
    end

end

function [f,fij,lam,mu] = A2_f_ods(n,d,u,a,gamma,m,f_last,fij_last,lam_last,mu_last,demand)
    
    lam = lam_last;
    f = f_last; 
    fij = fij_last;
    mu = mu_last;

    lam_prev = -10*ones(n);
    f_prev = -10*ones(n);
    fij_prev = -10*ones(n,n,n,n);
    mu_prev = -10*ones(n,n,n);

    beta_1 = 1e-2/gamma; %gamma 10
    beta_2 = 4e-1/gamma; %gamma 10

    max_dif = 1;
    cons_f_a  =20;
   

    while ( mean(mean(abs(f-f_prev)))/(beta_1*gamma) > 3e-2 ) ...
       || (cons_f_a  > 10) ...
       || (max_dif > 3e-1)  ...
          || ( mean(mean(max(max(abs(fij-fij_prev)))))/(beta_1*gamma) > 1e-1 )
        %|| (mean(mean((squeeze(sum(sum(fij,4),3)) < (a+1e-1)))) < 0.5) ...



        lam_prev = lam;
        f_prev =f;
        fij_prev = fij;
        mu_prev = mu;

       % disp(max_dif)
       %disp(['fij 1-2 = ',num2str(max(max(fij(:,:,1,2)))),', f 1-2 = ',num2str(f(1,2))]);
    

        mu = mu + beta_2.*(squeeze(sum( fij , 2 ))  - ...
          squeeze(permute(sum(fij , 1 ),[2 1 3 4])));

        for o=1:n
            for des=1:n
                mu(o,o,des) = mu(o,o,des) + beta_2.*(-f(o,des));
                mu_o(o,des) = mu(o,o,des);
                mu(des,o,des) = mu(des,o,des) + beta_2.*(f(o,des));
                mu_d(o,des) = mu(des,o,des);
            end
        end

      %  mu = max(mu,0);

        for o=1:n
            for des=1:n
                mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
                mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
            end
        end

        for i=1:n
            for j=1:n
                lam(i,j) = lam(i,j) + beta_2.*(   sum(sum(squeeze(permute(fij(i,j,:,:),[3,4,1,2]).*demand))) - a(i,j) );
            end
        end
        lam = max(lam,0);

        for i=1:n
            lam(i,i) = 100;
        end

        f = f - beta_1*(-m  + gamma.*demand.*(log(f./(1-f))  - u) - mu_o + mu_d );
        f = min(0.99,f);
        f = max(0.01,f);



        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( gamma*demand(o,des).*d + mu_i(:,:,o,des) - mu_j(:,:,o,des) + lam);
            end
        end

        fij = min(1,fij);
        fij = max(0,fij);

        cons_f_a = 0;
        for i=1:n
            for j=1:n
                if squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) >= (a(i,j) + 2e-1)
                    cons_f_a = cons_f_a + 1;
                end
            end
        end

        max_dif = 0;
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd))  );

                end
                
            end
        end
      %  disp(dif_ok./total_od);
       % disp(max_dif);

        for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
            mu(:,i,i) = 0;
        end

    end

end


function [f,fij] = A1_ods(n,d,u,a,f_last,fij_last,lam,mu)
    f = f_last;
    fij = fij_last; 

    mu_o = zeros(n);
    mu_d = zeros(n);

    mu_i = zeros(n,n,n,n);
    mu_j = mu_i;

    for o=1:n
        for des=1:n
            mu_o(o,des) = mu(o,o,des);
            mu_d(o,des) = mu(des,o,des);
        end
    end

    for o=1:n
        for des=1:n
            mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
            mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
        end
    end

    f_prev = -10*ones(n);
    fij_prev = -10*ones(n,n,n,n); %i,j,o,d

    beta_1 = 1e-2;

    q = 0;
    max_dif = 1;

    while ( max(max(abs(f-f_prev)))/(beta_1) > 4e-2 ) ...
   || ( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1) > 1e-2 ) ...
   || (max_dif > 7e-2)

        f_prev = f;
        fij_prev  = fij;

        % disp(max_dif);
        % disp( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1)   );

        f = f - beta_1*(log(f./(1-f))  - u - mu_o + mu_d );
        f = min(0.99,f);
        f = max(0.01,f);

        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( d + lam + mu_i(:,:,o,des) - mu_j(:,:,o,des) );
            end
        end
        fij = min(1,fij);
        fij = max(0,fij);

       for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
        end


        max_dif = 0;
        
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd)));
                end
                
            end
        end


    end
    disp('converge A1');


end


function [a_best,fg_best,fijg_best,lamg_best,mug_best,a_vec,f_vec,fij_vec,lamf_vec,muf_vec,obj_vec,obj_vec_f, his_time, time_best] = A3_ods(n,d,u,gamma,m,c,a_start,beta_base,demand)
    t_start = tic;
    a = a_start;

    fg = 0.01*ones(n);
    fijg = zeros(n,n,n,n);
    lamg = zeros(n);
    mug = zeros(n,n,n);

    ff = fg;
    fijf = fijg;
    lamf = lamg;
    muf = mug;

    a_prev = -10*ones(n);
    f_prev = -10*ones(n);
    fij_prev = -10*ones(n,n,n,n);
    lam_prev = -10*ones(n);
    mu_prev = -10*ones(n,n,n);

    % beta = 1e-4; %gamma 10
    % beta = 1e-5;
    % beta = 1e-6;
    
    beta = beta_base;

    numdatos = 3e2;
    obj_vec_lasts_g = zeros(1,numdatos);
    obj_vec_lasts_f = zeros(1,numdatos);
    
    hplot = plot(1:numdatos,obj_vec_lasts_g,'-b',1:numdatos,obj_vec_lasts_f,'--b','LineWidth',1.5);
    xlabel('Last 300 results');
    ylabel('F(a,f(a))');
    ylim([-10,inf]);
    %legend('F(a,f)','Obj(a,\theta)','Location','northwest');



    a_vec = [a];
    f_vec = [fg];
    fij_vec = [fijg];
    lamf_vec = [lamf];
   % lamg_vec = [lamg];
    muf_vec = [muf];
   % mug_vec = [mug];
    obj_vec = [];
    obj_vec_f = [];
    his_time = [];
    lamg_prev = -10*ones(n);
    max_iters = 2e5;
    obj_best = -1e6;
    k = 0;
    q = 0;
    control = 0;

    % while ((  max(max(abs(fg - f_prev))) > 1e-2) || ( sum(sum((abs(a-a_prev)./a_prev)/beta))  > 1) || ...
    %       (max(max(max(max(abs(fijg-fij_prev))))) > 1e-2) || ...
    %       (max(max(abs(lamf-lam_prev)))/beta > 1e-2) ||...
    %       (max(max(max(abs(muf-mu_prev))))/beta > 1e-2)) && ((k < max_iters)) %  ||  (obj <= (max(obj_vec) - 0.5 )) ) && (k < (max_iters + 2e4))

    while (k < max_iters)

        %disp(['a gap = ',num2str( sum(sum((abs(a-a_prev)./a_prev)/beta)) )]);

        
     %   disp( (abs(a-a_prev)./a_prev) );

     
        a_prev = a;
        f_prev = fg;
        ff_prev = ff;
        fij_prev = fijg;
        fijf_prev = fijf;
        lam_prev = lamf;
        lamg_prev = lamg;
        mu_prev = muf;
        mug_prev = mug;

        %disp(lamg);

        if (max(max(lamg_prev - 1000*eye(n))) < 1) && (k > 1)
            beta = 1e-2;
          %  disp('acelero');
        else
            % beta = 1e-4; %gamma 10
            % beta = 1e-5;
            % beta = 1e-6;
            beta = beta_base;
     
        end
       

        [fg,fijg,lamg,mug] = A2_ods(n,d,u,a,f_prev,fij_prev,lamg_prev,mug_prev,demand);
       % disp('converge el lower');
        [ff,fijf,lamf,muf] = A2_f_ods(n,d,u,a,gamma,m,ff_prev,fijf_prev,lam_prev,mu_prev,demand);
        gt = c + gamma.*lamg - lamf;
        a = a - beta.*gt;


         %  disp(['dlamg/dx = '   ,num2str(max(max(abs(lamg-lamg_prev)./abs(a-a_prev))))]);

        a = max(0.001,a);
      %  a = min(0.999,a);
      %  a_vec = [a_vec a(1,3)];
      %  f_vec = [f_vec fg];

      aa_ob = a;
      aa_ob(aa_ob < 1e-2) = 0;
      ff_ob = fg;
      ff_ob_f = ff;
      ff_ob(ff_ob < 2e-2) = 0;
      ff_ob_f(ff_ob_f < 2e-2) = 0;

        obj = sum(sum(m.*demand.*ff_ob)) - sum(sum(aa_ob.*c));
        obj_f = sum(sum(m.*demand.*ff_ob_f)) - sum(sum(aa_ob.*c));

        obj_vec = [obj_vec obj];
        obj_vec_f = [obj_vec_f obj_f];
            c_time=toc(t_start);
            his_time = [his_time, c_time];


        if (max(obj_vec) == obj)
            a_best = a;
            fg_best = fg;
            fijg_best = fijg;
            lamg_best = lamg;
            mug_best = mug;
            time_best = c_time;
        end

      %  lamf_vec = [lamf_vec lamf];
      %  lamg_vec = [lamg_vec lamg];
      %  muf_vec = [muf_vec muf];
     %   mug_vec = [mug_vec mug];

        % if (beta == 1e-2)
        %     k = k + 100;
        %     q = q + 100;
        %     control = control + 100;
        % else
            k = k+1;
            q = q+1;
            control = control + 1;
        % end
    %    if q >= 20000 gamma 10
         if (q >= 40000)  %gamma 20
%             disp('fijg = ');
%             disp(fijg);
%             disp('fijf = ');
%             disp(fijf);
            disp('fg = ');
            disp(fg);
            disp('a = ');
            disp(a);
            q = 0;


            % if ( max( obj_vec ) <= (obj_best + 0.01) )% && (k < max_iters))
            %      k = max_iters;
            % end
            obj_best = max(obj_vec);
            disp(['k = ',num2str(k),', best obj = ',num2str(obj_best)]);
        end

        if control >= 1e3
             disp(['k = ',num2str(k),', obj = ',num2str(obj),', beta = ',num2str(beta)]);
           %  ares = a;
           %  ares(ares < 0.1) = 0;
           %  fgres = fg;
           %  ffres = ff;
           %  fgres(fgres<0.02) = 0;
           %  ffres(ffres<0.02) = 0;
           %  res_g = sum(sum(m.*fgres.*demand))-sum(sum(c.*ares));
           %  res_f = sum(sum(m.*ffres.*demand))-sum(sum(c.*ares));
           %  obj_vec_lasts_g = [obj_vec_lasts_g(2:end),res_g];
           %  obj_vec_lasts_f = [obj_vec_lasts_f(2:end),res_f];
           % 
           %  set(hplot(1),'YData',obj_vec_lasts_g);
           %  set(hplot(2),'YData',obj_vec_lasts_f);
           %  drawnow;
           % % disp(a);
           % % disp(max(max(lamg-1000.*eye(n))));
            control = 0;
        end
        
       % disp(['k = ',num2str(k),', x = ',num2str(x),', y = ',num2str(yg), ', lam = ',num2str(lamf) , ', lam_op = ',num2str(lamg)]);
    end
   % disp(a);
   % disp(fg);
   % disp(fijg);
    %disp(['a = ',num2str(x),', y = ',num2str(yg), ', lam = ',num2str(lamf) , ', lam_op = ',num2str(lamg)]);
end



function [f,fij,lam,mu, error_f, error_fij] = ll_blocc(n,d,u,a,f_last,fij_last,lam_last,mu_last,demand, max_iter)
    lam = lam_last;
    f = f_last;
    fij = fij_last;
    mu = mu_last;

    mu_o = zeros(n);
    mu_d = zeros(n);

    mu_i = zeros(n,n,n,n);
    mu_j = mu_i;

    lam_prev = -10*ones(n); %i,j
    f_prev = -10*ones(n); %od
    fij_prev = -10*ones(n,n,n,n); %i,j,o,d
    mu_prev = -10*ones(n,n,n); %i,o,d

    beta_1 = 1e-2;
    beta_2 = 4e-1;
    q = 0;
    max_dif = 1;
    k=0;

    while ( mean(mean(abs(f-f_prev)))/(beta_1) > 2e-2 ) ...
       || ( mean(mean(mean(mean(abs(fij-fij_prev)))))/(beta_1) > 1e-2 ) ...
       || (cons_f_a > 10) ...
       || (max_dif > 2e-1) || (k<max_iter)
               % (mean(min((squeeze(sum(sum(fij,4),3)) < (a+3e-2)))) < 0.5) ... %sin demanda

        lam_prev = lam;
        f_prev = f;
        fij_prev  = fij;
        mu_prev = mu;
   
     %   disp(max_dif)

       

       mu = mu + beta_2.*(squeeze(sum( fij , 2 ))  - ...
          squeeze(permute(sum(fij , 1 ),[2 1 3 4])));

%         for o=1:n
%             for des=2
%                 mu(:,o,des) = mu(:,o,des) + beta*(sum( fij(:,:,o,des),2 ) - sum(fij(:,:,o,des) , 1)' );
%             end
%         end

        for o=1:n
            for des=1:n
                mu(o,o,des) = mu(o,o,des) + beta_2*(-f(o,des));
                mu_o(o,des) = mu(o,o,des);
                mu(des,o,des) = mu(des,o,des) + beta_2*(f(o,des));
                mu_d(o,des) = mu(des,o,des);
            end
        end

     %   mu = max(mu,0);

        for o=1:n
            for des=1:n
                mu_i(:,:,o,des) = mu(:,o,des)*ones(1,n);
                mu_j(:,:,o,des) = permute(mu_i(:,:,o,des),[2,1,3,4]);
            end
        end

        for i=1:n
            for j=1:n
                lam(i,j) = lam(i,j) + beta_2.*(   sum(sum(squeeze(permute(fij(i,j,:,:),[3,4,1,2]).*demand))) - a(i,j) );
            end
        end

%         for i=1:n
%             for j=1:n
%                 lam(i,j) = lam(i,j) + beta*(  sum(sum(squeeze(fij(i,j,:,:)))) - a(i,j));
%             end
%         end

        lam = max(0,lam);

        for i=1:n
            lam(i,i) = 100;
        end

        f = f - beta_1*(demand.*log(f./(1-f))  - demand.*u - mu_o + mu_d );

%         for o=1:n
%             for des=1:n
%                 f(o,des) = f(o,des) - beta*(log(f(o,des)./(1-f(o,des)))  - u(o,des) - mu(o,o,des) + mu(des,o,des) ); %puedo quitar algun bucle?
%             end
%         end

        
        f = min(0.99,f);
        f = max(0.01,f);

        for o=1:n
            for des=1:n
                fij(:,:,o,des) = fij(:,:,o,des) - beta_1*( demand(o,des).*d + lam + mu_i(:,:,o,des) - mu_j(:,:,o,des) );
            end
        end

%         for o=1:n
%             for des=1:n
%                 for i=1:n
%                     for j=1:n
%                         fij(i,j,o,des) = fij(i,j,o,des) - beta*(d(i,j) + lam(i,j) + mu(i,o,des) - mu(j,o,des)  );
%                     end
%                 end
%                % fij(:,:,o,des) = fij(:,:,o,des) - beta*(d + lam + mu_i - mu_j);
%             end
%         end

        fij = min(1,fij);
        fij = max(0,fij);

        for i=1:n
            f(i,i) = 0;
            fij(:,:,i,i) = 0;
            mu(:,i,i) = 0;
        end

        cons_f_a = 0;
        for i=1:n
            for j=1:n
                if squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) >= (a(i,j) + 2e-1)
                    cons_f_a = cons_f_a + 1;
                end
            end
        end



        max_dif = 0;
        
        for oo=1:n
            for dd=1:n
                if (oo ~= dd)
                    max_dif = max(max_dif, abs(sum(fij(oo,:,oo,dd),2) - f(oo,dd))  );
                end
                
            end
        end


        if q == 100
%             disp('fij 1,2 = ');
%             disp(fij(:,:,1,2));
%             disp('f 1,2 = ');
            %disp(f);
            q = 0;
        else
            q = q+1;
        end
        

        error_f=norm(f-f_last, 'fro');
        error_fij=norm(fij-fij_last, 'fro');
        disp(k);
        disp(error_f+error_fij);
        k=k+1;
        
        
    end

end



