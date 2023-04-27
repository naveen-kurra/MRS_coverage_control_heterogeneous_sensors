N = 5;
initial_conditions = [1.25, 0.25, 0;1, 0.5, 0;1, -0.5, 0;-1, -0.75, 0;0.1, 0.2, 0;0.2, -0.6, 0;-0.75, -0.1, 0;-1, 0, 0;-0.8, -0.25, 0;1.3, -0.4, 0];
r = Robotarium('NumberOfRobots', N, 'ShowFigure',true,InitialConditions = initial_conditions(1:N, :)');
x_min=-1.6;
x_max=1.6;
x_res=0.3;
y_min=-1;
y_max=1;
y_res=0.3;
c_value= 0.5;
epsilon = 0.5;
%%paramerters%%
sig = 2;
ann = 1;
psi = 2;

filename = 'data.csv';
final_goal_points = generate_initial_conditions(N, ...
    'Width', r.boundaries(2)-r.boundaries(1)-r.robot_diameter-0.25, ...
    'Height', r.boundaries(4)-r.boundaries(3)-r.robot_diameter-0.2, ...
    'Spacing', 0.5);
line_width = 5;
sum_C_x = 0;
sum_C_y = 0;
iterations = 1000;
wi_hi = zeros(iterations,N);
final_costs1 = zeros(iterations,1);
%     sensor_robots=[1 1 1 0 0;0 0 1 1 1];
%     Weights = [1 1 2 0 0;0 0 1 2 2];         %Scenario 4
%     Health = [1 1 1 0 0;0 0 1 1 1];

A = ones(N) - eye(N);
G = graph(A);
D = diag(sum(A));
L = D - A;
[~, uni_to_si_states] = create_si_to_uni_mapping();
si_to_uni_dyn = create_si_to_uni_dynamics_with_backwards_motion();
uni_barrier_cert_boundary = create_uni_barrier_certificate_with_boundary();


dxi = zeros(2, N);
CM = rand(N,3);

   for i=1:N
    robot_caption = sprintf('Robot %d', i);
    d(i) = plot(final_goal_points(1,i), final_goal_points(2,i),'s','LineWidth',line_width,'Color',CM(i,:));
    robot_labels{i} = text(500, 500, robot_caption, 'FontWeight', 'bold');
   end
fid = fopen(filename, 'w');
fprintf(fid, '%f,%f,%f,%f,%f,%f,%f\n', 'poses','weights','weigths-health','CW','distance_travelled', 'H', 'cum_dist');
cum_dist = zeros(iterations,1);
distances_travelled = zeros(iterations,N);

for t = 1:iterations
    costs_robots = zeros(N,1);
    x = r.get_poses();
    xi = uni_to_si_states(x);
    if t==1
        prev_location=xi;
    end
            for i=1:N
                c_v = zeros(N,2);
                Mass = zeros(N,1);
                sum_c_x = zeros(N,1);
                sum_c_y = zeros(N,1);
                all_sum_c_x = 0;
                all_sum_c_y = 0;
                    for qx = x_min:x_res:x_max
                        for qy = y_min:y_res:y_max
                               distances = zeros(N,1);
                                for robot = 1:N
                                    distances(robot) =  norm(([qx;qy]-xi(:,robot)) ,2);
                                end
                                [min_value, min_index] = min(distances);
                                c_v(min_index,1) = c_v(min_index,1) + (qx );
                                c_v(min_index,2) = c_v(min_index,2) + (qy);
                                Mass(min_index) = Mass(min_index) + 1 ;
                        end
                    end
                        if Mass(i) ~= 0
                            c_x = c_v(i,1)/Mass(i);
                            c_y = c_v(i,2)/Mass(i); 
                        else
                            c_x = 0;
                            c_y = 0; 
                        end
                        d(i).XData = c_x;
                        d(i).YData = c_y;
                        UiR = robustness(G,L,N,A,i,xi)
                        dxi(:, i) =ann*([c_x;c_y] - xi(:, i));
                        robot_labels{i}.Position = x(1:2, i);
            end
            prev_location= xi;
            norms = arrayfun(@(x) norm(dxi(:, x)), 1:N);
            threshold = 3/4*r.max_linear_velocity;
            to_thresh = norms > threshold;
            dxi(:, to_thresh) = threshold*dxi(:, to_thresh)./norms(to_thresh); 
            dxu = si_to_uni_dyn(dxi, x);
        %     dxu = uni_barrier_cert_boundary(dxu, x);
            r.set_velocities(1:N, dxu);   
            r.step();
end
fclose(fid);
r.debug();
function UiC = connectivity_maintenance(L,epsilon)
            pi = 0.5;
            [V, lambda2] = eigs(L, 2, 'sm');
            connectivity = lambda2(2,2);
            if connectivity > epsilon
                UiC = -(1-coth(connectivity - epsilon)^2)*1; % ∂(connectivity)/∂(pi)
            else
                UiC = 0;
            end
end
function UiR = robustness(G,L,N,A,i,xi)
            r = 0.5;
            m = numedges(G);
            edge_weights = zeros(m, 1);
            k = 1;
            for i = 1:numnodes(G)
                for j = i+1:numnodes(G)
                    if A(i,j) ~= 0
                        edge_weights(k) = L(i,i) + L(j,j) - 2*L(i,j);
                        k = k + 1;
                    end
                end
            end
            display(edge_weights)
            bc = centrality(G, 'betweenness', 'Cost',edge_weights);
            [gc, idx] = sort(bc, 'descend');
            
            j = 1;
            while ~isempty(conncomp(rmnode(G, idx(1:j)))) && j < numel(idx)
                j = j + 1;
            end
            if j == numel(idx)
                omega = 0;
            else
                omega = j;
            end
            Theta = omega/N;
%             display(Theta)
            hop1_neighbors = find(A(i,:));
            A2X = A * A;
            hop2_neighbors = A2X(i,:);
            hop2_neighbors = setdiff(find(hop2_neighbors), find(A(i,:)));
            num_hop1_neighbors = numel(hop1_neighbors);
            num_hop2_neighbors = numel(hop2_neighbors);
            prodb = num_hop1_neighbors + num_hop2_neighbors;
            beta = 1;
            A2 = A^2;
            neighbors2 = find(A2(i,:));
            display(neighbors2)
            reachable = unique([find(A(i,:)), neighbors2]);
            A_beta = A^beta;
            reachable_beta = find(any(A_beta(:,reachable),2));
            
            % Compute number of nodes at exactly 2 hops and reachable in at most beta hops
            nodes_2betahops = setdiff(reachable_beta, reachable);
            num_nodes = length(setdiff(reachable_beta, reachable));
            
            vulnerability_level = num_nodes/prodb;
            
            pos_sum_x = sum(xi(1,nodes_2betahops));
            pos_sum_x = pos_sum_x/num_nodes;
            pos_sum_y = sum(xi(2,nodes_2betahops));
            pos_sum_y = pos_sum_y/num_nodes;
            barycenter = [pos_sum_x;pos_sum_y];
            zzp = rand;
            if vulnerability_level > zzp
                eps = 1;
            else
                eps = 0;
            end
            disp(barycenter)
            UiR = eps*(barycenter - xi(:,i))/norm(barycenter-xi(:,i))*0.5;
            disp(UiR)
end 