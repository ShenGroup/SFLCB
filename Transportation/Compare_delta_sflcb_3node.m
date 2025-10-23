%% 


%% Inicializacion red 3 nodos

clear all;
close all;
clc;


n = 3;
d = 2.*(ones(n) - eye(n));
t = [0,1,10;1,0,2;10,2,0];
t_ext = 3.*(ones(n) - eye(n));
c = [0,1,10;1,0,3;10,3,0];
m = [0,2,6;2,0,1;6,1,0];
w = ones(n);

max_iter = 20000;
rho1 = 1000;
rho2 = 1000;
etax = 3e-4;
etay = 3e-4;
etaz = 3e-4;
etau = 3e-4;
etav = 3e-4;

delta_list = [0.01, 0.05, 0.1, 0.5, 1];
seeds = [1, 2, 3];
num_seeds = length(seeds);


results = struct();

for j = 1:length(delta_list)
    delta = delta_list(j);
    all_ob_z = []; 

    for s = 1:num_seeds
        rng(seeds(s));
        [~, ~, ~, his_ob_z, ~, ~, ~] = ...
            sflcb(n, c, m, t, t_ext, w, max_iter, delta, rho1, rho2, etax, etay, etaz, etau, etav);


        min_len = min(max_iter, length(his_ob_z));
        his_ob_z = his_ob_z(1:min_len);

        if size(all_ob_z, 1) < min_len
            all_ob_z(end+1:min_len, :) = NaN;
        end
        all_ob_z(:, end+1) = his_ob_z(:);
    end


    results(j).delta = delta;
    results(j).seeds = seeds;
    results(j).his_ob_z = all_ob_z;  % (iteration, seed)
end


save('his_ob_z_all_deltas.mat', 'results');




function [x_best, z_best, zij_best, his_ob_z, his_time, time_best, obj_best]=sflcb(n, c, m, t, t_ext, w, max_iter, delta, rho1, rho2, etax, etay, etaz, etau, etav)
    
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



        if all(hxq < 1e-2, 'all') && max(abs(gradmv_), [], 'all') < 2e-1
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
            % disp(mu_od)
    
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
