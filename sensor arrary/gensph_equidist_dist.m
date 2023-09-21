function [x,y,z,dist] = gensph_equidist_dist(r, n_points)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% generate points equidistantly spaced along the surpace of a sphere with
% radius "r", returns [x,y,z] cartesian points
% algorithm based on "How to generate equidistributed points on the surface
% of a sphere" by Markus Deserno
% modified to take in distance between points as an input

cartcor=[];
N_count=0;
area=(4*pi*r^2)/n_points;
d=sqrt(area);
%area=d^2;
M_theta=round(pi/d);
d_theta=pi/M_theta;
d_phi=area/d_theta;

for i=(1:M_theta)
    theta = pi*(i-1+0.5)/M_theta; %subtract 1 to adjust for matlab indexing
    M_phi=round(2*pi*sin(theta)/d_phi); 
    for j=(1:M_phi)
        phi=2*pi*(j-1)/M_phi;
        cartcor(N_count+1,:)=[r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)];
        N_count=N_count+1
    end
end

x=cartcor(:,1);
y=cartcor(:,2);
z=cartcor(:,3);
dist=sqrt(abs((x(2)-x(1)^2)+(y(2)^2-y(1)^2)+(z(2)^2-z(1)^2)));
end