

%% Inicializacion red 3 nodos

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

 t = distance;

 t_ext = [0,1.6,0.8,2,1.6,2.5,3,2.5,0.8; ...
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

w = 1e-1.*[0,9,26,19,13,12,13,8,11;
          11,0,14,26,7,18,3,6,12;
          30,19,0,30,24,8,15,12,5;
          21,9,11,0,22,16,25,21,23;
          14,14,8,9,0,20,16,22,21;
          26,1,22,24,13,0,16,14,12;
          8,6,9,23,6,13,0,11,11;
          9,2,14,20,18,16,11,0,4;
          8,7,11,22,27,17,8,12,0];



max_iter = 300000;
delta=0.25;
rho1=50;
rho2=50;
etax=3e-5;
etay=3e-5;
etaz=3e-5;
etau=3e-5;
etav=3e-5;


rng(1);




[x_best, z_best, zij_best, his_ob_y, his_ob_z, his_ob_diff, his_time, time_best, obj_best] = sflcb(n, c, m, t, t_ext, w, max_iter, delta, rho1, rho2, etax, etay, etaz, etau, etav);

plot(his_time, his_ob_z, '-o');
xlabel('Time (s)');
ylabel('UL utility');
title('9 Node');
grid on;


function [x_best, z_best, zij_best, his_ob_y, his_ob_z, his_ob_diff, his_time, time_best, obj_best]=sflcb(n, c, m, t, t_ext, w, max_iter, delta, rho1, rho2, etax, etay, etaz, etau, etav)
    
    t_start = tic;
    his_ob_diff = zeros(max_iter, 1);
    his_ob_y = zeros(max_iter, 1);
    his_ob_z = zeros(max_iter, 1);
    his_time = zeros(max_iter, 1);
    obj_best = -1e6;

    x = 0.7.*rand(n,n) + 0.3;

    y = 0.01*ones(n); z = y;
    p = zeros(n,n,n,n); q = p;

    mv = zeros(n,n,n); mu = mv;

    u = zeros(n,n); v = u;

    y(1:n+1:end) = 0;
    z(1:n+1:end) = 0;

    alpha = compute_h(x,p,w); beta = alpha;


    for iter = 1:max_iter
        % compute gradients
        hxp = compute_h(x, p, w);
        hxq = compute_h(x, q, w);
        h_alpha=hxp - alpha;
        h_beta=hxq - beta;

        gradx = delta * c - u + v - rho1 * h_alpha + rho2 * h_beta;



        gradalpha = - u - rho1 * h_alpha;
        gradbeta = v + rho2 * h_beta;




        gradu = h_alpha;
        gradv = -h_beta;

        gradmu = compute_l(y, p, n);
        gradmv_ = compute_l(z, q, n);




        grady = grad_y(y, m, mu, rho1, w, t_ext, delta, gradmu);
        gradz = grad_z(z, mv, rho2, w, t_ext, gradmv_);
        gradp = grad_p(p, mu, u, rho1, w, t, gradmu, h_alpha);
        gradq = grad_q(q, mv, v, rho2, w, t, gradmv_, h_beta);



        x = x - etax * gradx;
        y = y - etay * grady;
        z = z + etaz * gradz;
        p = p - etay * gradp;
        q = q + etaz * gradq;
        alpha = alpha - etay * gradalpha;
        beta = beta + etaz * gradbeta;

        u = u + etau * gradu;
        v = v - etav * gradv;
        mu = mu + etau * gradmu;
        mv = mv + etav * gradmv_;


        x = max(0.001,x);

        y = min(0.99,y);
        y = max(0.01,y);
        z = min(0.99,z);
        z = max(0.01,z);

        p = min(1,p);
        p = max(0,p);
        q = min(1,q);
        q = max(0,q);

        alpha=min(0,alpha);
        beta=min(0,beta);

        mu=max(0,mu);
        mv=max(0,mv);

        % Set diagonals of y, z to zero
        y(1:n+1:end) = 0;
        z(1:n+1:end) = 0;

        
        for i = 1:n
            p(i,i,:,:) = 0;
            q(i,i,:,:) = 0;
            mu(:, i, i) = 0;
            mv(:, i, i) = 0;
        end



        obj_y = sum(sum(m.*y)) - sum(sum(x.*c));
        obj_z = sum(sum(m.*z)) - sum(sum(x.*c));

        his_ob_y(iter)=obj_y;
        his_ob_z(iter)=obj_z;
        his_ob_diff(iter)=abs(obj_z-obj_y);
        his_time(iter)=toc(t_start);

        count= sum(hxq < 0.2, 'all');

        if count >= n^2 - 10 && max(abs(gradmv_), [], 'all') < 2e-1
            disp([iter, obj_z]);
            if obj_z > obj_best
                disp(['best ', num2str(obj_z)]);
                obj_best = obj_z;
                x_best = x;
                z_best = z;
                zij_best = q;
                time_best = his_time(iter);
            end

        end
    end
end



function h = compute_h(x, p, w)
    h = squeeze(sum(sum(w .* p, 1), 2)) - x;
end


function l = compute_l(y, p, n)

    l = squeeze(sum(p, 4)) - squeeze(sum(p, 3));  % size: (n, n, n)


    for o = 1:n
        for d = 1:n
            l(o,d,o) = l(o,d,o) - y(o,d);
            l(o,d,d) = l(o,d,d) + y(o,d);
        end
    end
end



function grad = grad_y(y, m, mu, rho1, w, t_ext, delta, l)
    n = size(y,1);
    eps_ = 1e-6;
    grad = -delta * m + w .* (log(y + eps_) - log(1 - y + eps_)) - w .* t_ext;
    
    %l = compute_l(y, p, n);
    for o = 1:n
        for d = 1:n
            grad(o,d) = grad(o,d) ...
                        - mu(o,d, o) - rho1 * l(o,d, o) ...
                        + mu(o,d, d) + rho1 * l(o,d, d);
        end
    end
end

function grad = grad_z(z, mv, rho2, w, t_ext, l)
    n = size(z,1);
    eps_ = 1e-6;
    grad = - w .* (log(z + eps_) - log(1 - z + eps_)) + w .* t_ext;
    
    %l = compute_l(z, q, n);
    for o = 1:n
        for d = 1:n
            grad(o,d) = grad(o,d) ...
                        + mv(o,d, o) + rho2 * l(o,d, o) ...
                        - mv(o,d, d) - rho2 * l(o,d, d);
        end
    end
end




function grad = grad_p(p, mu, u, rho1, w, t, l, h_alpha)
    n = size(p,1);

    grad = zeros(n,n,n,n);
    base_term = t + u + rho1 * h_alpha;  % size: (n,n)
    
    for o = 1:n
        for d = 1:n
            w_od = w(o,d);
            mu_od = squeeze(mu(o,d, :));  % size: (n,1)
            l_od  = squeeze(l(o,d, :));   % size: (n,1)

    
            % Broadcasting over (i,j)
            grad(o,d,:,:) = w_od * base_term + ...
                            mu_od - mu_od' + ...
                            rho1 * (l_od - l_od');
        end
    end
end

function grad = grad_q(q, mv, v, rho2, w, t, l, h_beta)
    n = size(q,1);

    grad = zeros(n,n,n,n);
    base_term = -t - v - rho2 * h_beta;  % size: (n,n)
    
    for o = 1:n
        for d = 1:n
            w_od = w(o,d);
            mv_od = squeeze(mv(o,d, :));  % size: (n,1)
            l_od  = squeeze(l(o,d, :));   % size: (n,1)
    
            % Broadcasting over (i,j)
            grad(o,d,:,:) = w_od * base_term + ...
                            -mv_od + mv_od' - ...
                            rho2 * (l_od - l_od');
        end
    end
end
