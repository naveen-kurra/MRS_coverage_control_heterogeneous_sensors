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
final_costs = zeros(iterations,1);

    sensor_robots=[1 1 1 0 0;0 0 1 1 1];
    Weights = [1 1 2 0 0;0 0 1 2 2];         %Scenario 4
    Health = [1 1 1 0 0;0 0 1 1 1];
% 
%     sensor_robots=[1 1 1 0 0;0 0 1 1 1];
%     Weights = [1 1 1 0 0;0 0 1 1 1];        %Scenario 5
%     Health = [1 1 2 0 0;0 0 1 2 1];

%     sensor_robots=[1 1 0 0 0 0 0 0 0 0;0 0 1 1 0 0 0 0 0 0;0 0 0 0 1 1 0 0 0 0;0 0 0 0 0 0 1 1 0 0;0 0 0 0 0 0 0 0 1 1];
%     Weights = [1 2 0 0 0 0 0 0 0 0;0 0 2 1 0 0 0 0 0 0;0 0 0 0 1 2 0 0 0 0;0 0 0 0 0 0 2 1 0 0;0 0 0 0 0 0 0 0 1 2];         %Scenario 3
%     Health = [1 2 0 0 0 0 0 0 0 0;0 0 2 1 0 0 0 0 0 0;0 0 0 0 1 2 0 0 0 0;0 0 0 0 0 0 2 1 0 0;0 0 0 0 0 0 0 0 1 2];

%     sensor_robots=[1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];
%     Weights = [1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];   %Scenario 2
%     Health = [1 1 2 1 1 0 0 0 0 0;0 0 0 0 0 2 1 1 1 2];

%     sensor_robots=[1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];
%    Weights = [1 1 1 1 1 0 0 0 0 0;0 0 0 0 0 1 1 1 1 1];        %Scenario 1
%     Health = [1 1 1 1 1 1 1 1 1 1];


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
fid = fopen(filename, 'w');
fprintf(fid, '%f,%f,%f,%f,%f,%f,%f\n', 'poses','weights','weigths-health','CW','distance_travelled', 'H', 'cum_dist');
cum_dist = zeros(iterations,1);
distances_travelled = zeros(iterations,N);
dist_robot = zeros(iterations,1);
% check1=zeros(iterations,N,size(sensor_robots,1));
w_h = zeros(iterations,2);
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
                                distances2 = zeros(N,1);
                                for robot = 1:length(index)
%                                     distances(index(robot)) =  norm(([qx;qy]-xi(:,index(robot))) ,2)-Weights(j,index(robot));
                                    distances(index(robot)) =  (norm(([qx;qy]-xi(:,index(robot))) ,2));
%                                     display(Weights(j,index(robot)));
                                    distances2(index(robot)) =  (norm(([qx;qy]-xi(:,index(robot))) ,2))^2-Weights(j,index(robot));
                                end
                                [min_value, min_index] = min(distances(index));
                                [min_value2, min_index2] = min(distances2(index));
                                min_index = find(distances == min_value);
                                min_index2 = find(distances2 == min_value2);
                                c_v(min_index,1) = c_v(min_index,1) + (qx * get_sensor);
                                c_v(min_index,2) = c_v(min_index,2) + (qy * get_sensor);
                                Mass(min_index) = Mass(min_index) + get_sensor ;
                                final_costs(t)=final_costs(t)+0.5*((min_value2))*get_sensor;    
%                                 display(Mass(min_index));
                        end
                    end
                    if sensor_robots(j,i)==1
                      
                       
                        if Mass(i) ~= 0
                            
                            c_x = c_v(i,1)/Mass(i);
                            c_y = c_v(i,2)/Mass(i); 
                        else
                            c_x = 0;
                            c_y = 0;
                            
                        end
                        
                        dist = norm(([c_x;c_y]-xi(:,i)) ,2);
%                         costs_robots(index(p))=costs_robots(index(p))+0.5*((dist)^2 - Weights(j,index(p)))*get_sensor;
                        d(i).XData = c_x;
                        d(i).YData = c_y; 
                        dxi(:, i) = ([c_x;c_y] - xi(:, i));

%                         

                        

%                         display(Weights(j,i));
                       
                    end
                    
%                      for p=1:length(index) 
%                      if sensor_robots(j,i) == 1
%                         
%                         
%                      end
                     
                     distances_travelled(t,i) = distances_travelled(t,i)+norm((prev_location(:,i)-xi(:,i)) ,2);
                     cum_dist(t) = cum_dist(t) + distances_travelled(t,i);
                    
                    
                    robot_labels{i}.Position = x(1:2, i);
            end
                     
        end
        prev_location = xi;
%         for m=1:N
%             final_costs(t) = final_costs(t)+costs_robots(m);
%         end
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
            fprintf(fid, '%f,%f,%f,%f,%f,%f,%f\n', xi,[c_x,c_y],distances_travelled, final_costs(t), cum_dist(t));
end
fclose(fid);
r.debug();

function sensor_value = get_sensor()
            sensor_value = 1;
end