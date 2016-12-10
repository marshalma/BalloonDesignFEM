function ShellDynamic(obj_file_path)

global video;
video = [];
E = 0.1;
nu = 0.45;
mu = E * nu / ((1+nu)*(1-2*nu));
lambda = 0.5 * E / (1 + nu);
rho = 1.0;
pressure = 0;
damping = 2;
% h = 1e-5;


% Creates triangles from a regular grid of nodes
[nodes,tris] = read_vertices_and_faces_from_obj_file(obj_file_path);
nNodes = length(nodes);
nTris = length(tris);
draw(0, nodes, tris);

% Compute parameters for each element
elements = [];
mass = zeros(nNodes, 1);
for k = 1 : nTris
	e12 = nodes(tris(k,2),:) - nodes(tris(k,1),:);
    u = e12 / norm(e12, 2);
    e13 = nodes(tris(k,3),:) - nodes(tris(k,1),:);
    d = cross(u,e13);
    d = d / norm(d, 2);
    v = cross(d,u);
    
    elements(k).e12 = e12;
    elements(k).e13 = e13;
    elements(k).u = u;
    elements(k).v = v;
    elements(k).d = d;
    elements(k).T = [u' v' d'];
    elements(k).A = 0.5 * norm(cross(e12,e13),2);
    
    mass(tris(k,1)) = mass(tris(k,1))+ elements(k).A * rho / 3;
    mass(tris(k,2)) = mass(tris(k,2))+ elements(k).A * rho / 3;
    mass(tris(k,3)) = mass(tris(k,3))+ elements(k).A * rho / 3;
end

% if ~isempty(video)
% 	video.close();
% end

t0 = -inf;
drawHz = 100;
dt = 0.01;
tEnd = 3;
velo = zeros(nNodes, 3);
for t = 0 : dt : tEnd
    t
	% Draw scene
    if t - t0 > 1 / drawHz
		draw(t,nodes,tris);
		t0 = t;
    end
    
    pforce = zeros(nNodes, 3);
    ppforce = zeros(nNodes, 3);
%     for i = 1:nTris
%         % calculating pressure force for each node by iterating all element
%         e12 = nodes(tris(i,2),:) - nodes(tris(i,1),:);
%         e13 = nodes(tris(i,3),:) - nodes(tris(i,1),:);
%         d = cross(e12, e13);
%         d = d / norm(d,2);
%         area = 0.5 * norm(cross(e12,e13),2);
%         force = pressure * d * area;
%         pforce(tris(i,1),:) = pforce(tris(i,1),:) + force / 3;
%         pforce(tris(i,2),:) = pforce(tris(i,2),:) + force / 3;
%         pforce(tris(i,3),:) = pforce(tris(i,3),:) + force / 3;
%     end
    
    for i = 1:nTris
        % calculating elastic force for each node by interating all
        % element
        e_12 = nodes(tris(i,2),:) - nodes(tris(i,1),:);
        e_13 = nodes(tris(i,1),:) - nodes(tris(i,3),:);
        e_23 = nodes(tris(i,3),:) - nodes(tris(i,2),:);
        dd = cross(e_12, e_13);
        dd = dd / norm(dd,2);
        norm12 = cross(dd, e_12);
        norm23 = cross(dd, e_23);
        norm13 = cross(dd, e_13);
        
        F = [e_12' e_13' dd'] * inv([elements(i).e12' -elements(i).e13' elements(i).d']) * elements(i).T;
        epsilon = 0.5 * (F' * F - eye(3));
        P = F * (2*mu*epsilon + lambda*trace(epsilon)*eye(3));
        sigma = P * F / det(F);
        f12 = sigma * norm12';
        f23 = sigma * norm23';
        f13 = sigma * norm13';
%         ppforce(tris(i,1),:) = ppforce(tris(i,1),:) + f12' / 2 + f13' / 2;
%         ppforce(tris(i,2),:) = ppforce(tris(i,2),:) + f12' / 2 + f23' / 2;
%         ppforce(tris(i,3),:) = ppforce(tris(i,3),:) + f23' / 2 + f13' / 2;
        pforce(tris(i,1),:) = pforce(tris(i,1),:) + f12' / 2 + f13' / 2;
        pforce(tris(i,2),:) = pforce(tris(i,2),:) + f12' / 2 + f23' / 2;
        pforce(tris(i,3),:) = pforce(tris(i,3),:) + f23' / 2 + f13' / 2;
    end
	
	% Integrate velocity and position
	% ### TODO ###
    quiver3(nodes(:,1),nodes(:,2),nodes(:,3),ppforce(:,1),ppforce(:,2),ppforce(:,3));
    for k = 1 : nNodes
        % using implicit damping
        velo(k,:) = (mass(k) * velo(k,:) + dt * pforce(k,:)) / (mass(k) + dt * damping * mass(k));
        % nodes(k).v = (nodes(k).m * nodes(k).v + dt * nodes(k).f) / (nodes(k).m + dt * damping * nodes(k).m);
        nodes(k,:) = nodes(k,:) + velo(k,:) * dt;
    end
end

if ~isempty(video)
	video.close();
end


end





%%
function draw(t,nodes,tris)

global video;

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
    zlabel('Z');
	axis equal;
%   	axis(3*[-1 1 -1 1 -1 1]); % Change axis limits here
	grid on;
	view(3);
	colormap jet;
	caxis([0 25]); % Change color limits here (comment out for auto)
	cb = colorbar;
	ylabel(cb, 'stress')
	video = VideoWriter('output','MPEG-4');
	video.open();
end
cla;
hold on;

figure(1);
trimesh(tris, nodes(:,1), nodes(:,2), nodes(:,3));
str = sprintf('t = %.4f', t);
title(str);
rotate3d on;



% quiver3(center_points(:,1),center_points(:,2),center_points(:,3),normal_vec(:,1),normal_vec(:,2),normal_vec(:,3));
% rotate3d on;
frame = getframe(gcf);
video.writeVideo(frame);

end


