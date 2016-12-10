function ShellStatic(obj_file_path)

% global video;
% video = [];
E = 100;
nu = 0.4;
global mu
mu = E * nu / ((1+nu)*(1-2*nu));
global lambda
lambda = 0.5 * E / (1 + nu);
global pressure
pressure = 0;
global h
h = 1e-5;


% Creates triangles from a regular grid of nodes
[nodes,tris] = read_vertices_and_faces_from_obj_file(obj_file_path);
nNodes = length(nodes);
nTris = length(tris);
draw(0, nodes, tris);

% Compute parameters for each element
elements = [];
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
    elements(k).A = 0.5 norm(cross(e12,e13),2);
end

% if ~isempty(video)
% 	video.close();
% end

opt = optimoptions('lsqnonlin',...
            'Jacobian','off',...
            'DerivativeCheck','off',...
            'Display','off');
lb = -inf(size(nodes));
ub =  inf(size(nodes));
[nodes_new,residue] = lsqnonlin(@(nodes_new)objFun(nodes_new, tris, elements), nodes, lb, ub, opt);
residue
draw(0, nodes_new, tris);

end


function pforce = objFun(nodes_new, tris, elements)
%     draw(0, nodes_new, tris);
    global pressure
    global mu
    global lambda
    global h
    pforce = zeros(length(tris), 3);
    for i = 1:length(tris)
        e12 = nodes_new(tris(i,2),:) - nodes_new(tris(i,1),:);
        e13 = nodes_new(tris(i,3),:) - nodes_new(tris(i,1),:);
        d = cross(e12, e13);
        d = d / norm(d,2);
        area = 0.5 * norm(cross(e12,e13),2);
        pforce(tris(i,1),:) = pforce(tris(i,1),:) + pressure * d * area / 3;
        pforce(tris(i,2),:) = pforce(tris(i,2),:) + pressure * d * area / 3;
        pforce(tris(i,3),:) = pforce(tris(i,3),:) + pressure * d * area / 3;
        
        e12 = nodes_new(tris(i,2),:) - nodes_new(tris(i,1),:);
        e13 = nodes_new(tris(i,3),:) - nodes_new(tris(i,1),:);
        e23 = nodes_new(tris(i,3),:) - nodes_new(tris(i,2),:);
        d = cross(e12, e13);
        d = d / norm(d,2);
        norm12 = cross(d, e12);
        norm23 = cross(d, e23);
        norm13 = -cross(d, e13);
        norm12 = norm12 / norm(norm12,2);
        norm23 = norm23 / norm(norm23,2);
        norm13 = norm13 / norm(norm13,2);
        
        F = [e12' e13' d'] * inv([elements(i).e12' elements(i).e13' elements(i).d']) * elements(i).T;
        epsilon = 0.5 * (F' * F - eye(3));
        P = F * (2*mu*epsilon + lambda*trace(epsilon)*eye(3));
        sigma = P * F / det(F);
        f12 = sigma * norm12';
        f23 = sigma * norm23';
        f13 = sigma * norm13';
        pforce(tris(i,1),:) = pforce(tris(i,1),:) + f12' / 2 + f13' / 2;
        pforce(tris(i,2),:) = pforce(tris(i,2),:) + f12' / 2 + f23' / 2;
        pforce(tris(i,3),:) = pforce(tris(i,3),:) + f23' / 2 + f13' / 2;
    end
%     e0 = elasticEnergy(nodes_new, tris, elements);
%     for i = 1:length(nodes)
%         diff = 0.01;
%         nodes_x_diff = nodes_new;
%         nodes_x_diff(i,1) = nodes_x_diff(i,1) + diff;
%         e_diff_x = elasticEnergy(nodes_x_diff, tris, elements);
%         nodes_y_diff = nodes_new;
%         nodes_y_diff(i,2) = nodes_y_diff(i,2) + diff;
%         e_diff_y = elasticEnergy(nodes_y_diff, tris, elements);
%         nodes_z_diff = nodes_new;
%         nodes_z_diff(i,3) = nodes_z_diff(i,3) + diff;
%         e_diff_z = elasticEnergy(nodes_z_diff, tris, elements);
%         gradient = ([e_diff_x, e_diff_y, e_diff_z] - e0) / diff;
%         res(i) = norm(pforce(i,:)-gradient, 2);
%     end
end

% function res = elasticEnergy(nodes_new, tris, elements)
%     
%     res = 0;
%     for i = 1:length(tris)
%         e12 = nodes_new(tris(i,2),:) - nodes_new(tris(i,1),:);
%         e13 = nodes_new(tris(i,3),:) - nodes_new(tris(i,1),:);
%         d = cross(e12, e13);
%         d = d / norm(d,2);
%         F = [e12' e13' d'] * inv([elements(i).e12' elements(i).e13' elements(i).d']) * elements(i).T;
%         epsilon = 0.5 * (F' * F - eye(3));
%         P = F * (2*mu*epsilon + lambda*trace(epsilon)*eye(3));
%         sigma = P * F / det(F);
%         res = res + sigma * h * elements(i).A;
%     end
% end

%%
function draw(t,nodes,tris)

% global video;

if t == 0
	clf;
	xlabel('X');
	ylabel('Y');
    zlabel('Z');
	axis equal;
%  	axis(1.5*[-1 1 -1 1 -1 1]); % Change axis limits here
	grid on;
	view(3);
	colormap jet;
	caxis([0 25]); % Change color limits here (comment out for auto)
	cb = colorbar;
	ylabel(cb, 'stress')
% 	video = VideoWriter('output','MPEG-4');
% 	video.open();
end
cla;
hold on;

figure(1);
trimesh(tris, nodes(:,1), nodes(:,2), nodes(:,3));
rotate3d on;

% center_points = zeros(length(tris), 3);
% normal_vec = zeros(length(tris), 3);
% for i = 1:length(tris)
%     e12 = nodes(tris(i,2),:) - nodes(tris(i,1),:);
%     e13 = nodes(tris(i,3),:) - nodes(tris(i,1),:);
%     d = cross(e12, e13);
%     d = d / norm(d,2);
%     normal_vec(i,:) = d;
%     center_points(i,:) = (nodes(tris(i,1),:) + nodes(tris(i,2),:) + nodes(tris(i,3),:)) / 3;
% end
% 
% hold on
% quiver3(center_points(:,1),center_points(:,2),center_points(:,3),normal_vec(:,1),normal_vec(:,2),normal_vec(:,3));
% rotate3d on;
% frame = getframe(gcf);
% video.writeVideo(frame);

end


