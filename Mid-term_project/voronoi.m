% Compute Voronoi diagram of the points
points = rand(5,2);
[V,C] = voronoin(points);

% Plot Voronoi diagram
figure
for i = 1:length(C)
    if all(C{i}~=1)   % If at least one index is 1, then it is an open region and we can't patch that
        patch(V(C{i},1),V(C{i},2),i); % use color i.
    end
end
plot(points(:,1),points(:,2),'k.','MarkerSize',15)
axis([0 1 0 1])
axis equal