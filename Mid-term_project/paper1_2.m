N = 5;  
initial_conditions = [1.25, 0.25, 0;1, 0.5, 0;1, -0.5, 0;-1, -0.75, 0;0.1, 0.2, 0;0.2, -0.6, 0;-0.75, -0.1, 0;-1, 0, 0;-0.8, -0.25, 0;1.3, -0.4, 0];
r = Robotarium('NumberOfRobots', N, 'ShowFigure',true,InitialConditions = initial_conditions(1:N, :)');
x_min=-1.6;
x_max=1.6;
x_res=0.15;
y_min=-1;
y_max=1;
y_res=0.15;
c_value= 0.5;
filename = 'dataS1.csv';
final_goal_points = generate_initial_conditions(N, ...
    'Width', r.boundaries(2)-r.boundaries(1)-r.robot_diameter-0.25, ...
    'Height', r.boundaries(4)-r.boundaries(3)-r.robot_diameter-0.2, ...
    'Spacing', 0.5);
line_width = 5;
sum_C_x = 0;
sum_C_y = 0;
iterations = 1000;
k=1;
wi_hi = zeros(iterations,N);
final_costs3 = zeros(iterations,1);
    sensor_robots=[1 1 1 0 0;0 0 1 1 1];
    Weights = [1 1 2 0 0;0 0 1 2 2];         %Scenario 4
    Health = [1 1 1 0 0;0 0 1 1 1];
% % 
%     sensor_robots=[1 1 1 0 0;0 0 1 1 1];
%     Weights = [1 1 1 0 0;0 0 1 1 1];         %Scenario 5
%     Health = [1 1 2 0 0;0 0 1 2 1];
% 
%     sensor_robots=[1 1 0 0 0 0 0 0 0 0;0 0 1 1 0 0 0 0 0 0;0 0 0 0 1 1 0 0 0 0;0 0 0 0 0 0 1 1 0 0;0 0 0 0 0 0 0 0 1 1];
%     Weights = [1 1 0 0 0 0 0 0 0 0;0 0 1 1 0 0 0 0 0 0;0 0 0 0 1 1 0 0 0 0;0 0 0 0 0 0 1 1 0 0;0 0 0 0 0 0 0 0 1 1];         %Scenario 3
%     Health = [1 2 0 0 0 0 0 0 0 0;0 0 2 1 0 0 0 0 0 0;0 0 0 0 1 2 0 0 0 0;0 0 0 0 0 0 2 1 0 0;0 0 0 0 0 0 0 0 1 2];
% %     Weights = [1 2 0 0 0 0 0 0 0 0;0 0 2 1 0 0 0 0 0 0;0 0 0 0 1 2 0 0 0 0;0 0 0 0 0 0 2 1 0 0;0 0 0 0 0 0 0 0 1 2];  

%     sensor_robots=[1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];
%     Weights = [1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];         %Scenario 2
%     Health = [1 1 2 1 1 0 0 0 0 0;0 0 0 0 0 2 1 1 1 2];

%     sensor_robots=[1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];
%     Weights = [1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];         %Scenario 1
%     Health = [1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];

filenames = {'w_h.csv','centeroids.csv','distance_tarvelled.csv','Hcosts.csv','cum_dist.csv'};
L = completeGL(N);
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
filepath = 'C:/Users/kurra/OneDrive/Desktop/multi_robot__systems/robotarium-matlab-simulator/Scenario1/';
cum_dist = zeros(iterations,1); 
distances_travelled = zeros(iterations,N);
dist_robot = zeros(iterations,1);
check1=zeros(iterations,N,size(sensor_robots,1));
w_h = zeros(iterations,2);
final_cent=zeros(iterations,N,2);
for t = 1:iterations
    costs_robots = zeros(N,1);
    x = r.get_poses();
    xi = uni_to_si_states(x);
    
    if t<2   
        prev_location = xi;
    end
        for j = 1:size(sensor_robots,1) % S is an array of sensors
            index=find(sensor_robots(j,:)==1);
            temp_weights = zeros(N,2);
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
                                for robot = 1:length(index)
%                                     distances(index(robot)) =  norm(([qx;qy]-xi(:,index(robot))) ,2)-Weights(j,index(robot));
                                    distances(index(robot)) = norm(([qx;qy]-xi(:,index(robot))) ,2)^2-Weights(j,index(robot));
%                                     display(Weights(j,index(robot)));
                                end
                                [min_value, min_index] = min(distances(index));
                                min_index = find(distances == min_value);
                                c_v(min_index,1) = c_v(min_index,1) + (qx * get_sensor);
                                c_v(min_index,2) = c_v(min_index,2) + (qy * get_sensor);
                                Mass(min_index) = Mass(min_index) + get_sensor ;
                                final_costs3(t)=final_costs3(t)+0.5*((min_value))*get_sensor;
%                                 display(Mass(min_index));
                        end
                    end
                    if sensor_robots(j,i)==1
                      neighbors = topological_neighbors(L, i);
                       temp_weight=0;
                        for f = neighbors
                            
                            temp_weight = temp_weight - ((Weights(j,i) - Weights(j,f)) - (Health(j,i) - Health(j,f)));
                            
                        end
                        if Mass(i) ~=0
                            temp_weight = 0.5*(k/Mass(i))*temp_weight*1;
                        else
                            temp_weight = 0; 
                        end
                        temp_weights(j,i) = temp_weight;
                       
                        if Mass(i) ~= 0
                            
                            c_x = c_v(i,1)/Mass(i);
                            c_y = c_v(i,2)/Mass(i); 
                        else
                            c_x = 0;
                            c_y = 0;
                            
                        end
                        
                        final_cent(t,i,1)=c_x;
                        final_cent(t,i,2)=c_y;

                        dist = norm(([c_x;c_y]-xi(:,i)) ,2);
%                         costs_robots(index(p))=costs_robots(index(p))+0.5*((dist)^2 - Weights(j,index(p)))*get_sensor;
                        d(i).XData = c_x;
                        d(i).YData = c_y; 
                        dxi(:, i) = ([c_x;c_y] - xi(:, i));

%                         display(Weights(j,i))
                        w_h(t,i) = Weights(j,i)-Health(j,i);
                        check1(t,i,j)=Weights(j,i); %Weights(j,i);

                        

%                         display(Weights(j,i));
                       
                    end
                     distances_travelled(t,i) = distances_travelled(t,i)+norm((prev_location(:,i)-xi(:,i)) ,2);
                     cum_dist(t) = cum_dist(t) + distances_travelled(t,i);
                    
                    
                    robot_labels{i}.Position = x(1:2, i);
            end
            for p=1:length(index)
                       Weights(j,index(p)) = Weights(j,index(p)) + temp_weights(j,index(p));

                    end
                     
        end
        
        

        if t>1
            dist_robot(t)=dist_robot(t-1)+cum_dist(t);
        end
            norms = arrayfun(@(x) norm(dxi(:, x)), 1:N);
            threshold = 3/4*r.max_linear_velocity;
            to_thresh = norms > threshold;
            dxi(:, to_thresh) = threshold*dxi(:, to_thresh)./norms(to_thresh); 
            dxu = si_to_uni_dyn(dxi, x);
        %     dxu = uni_barrier_cert_boundary(dxu, x);
            r.set_velocities(1:N, dxu);   
            r.step();
            datacent=[final_cent(t,:,1),final_cent(t,:,2)];
            variables = {w_h(t), datacent, distances_travelled(t,:),final_costs3(t),cum_dist(t)};
            for i = 1:length(variables)
                 dlmwrite(char(fullfile(filepath, filenames(i))),variables(i), '-append', 'delimiter', ',');
            end
%             csvwrite(fullfile(filepath, filename), data);
%             fprintf(fid, '%f,%f,%f,%f,%f,%f\n', xi,w_h(t),final_cent,distances_travelled(t), final_costs3(t),cum_dist(t));
end
r.debug();

function sensor_value = get_sensor()
            sensor_value = 1;
end