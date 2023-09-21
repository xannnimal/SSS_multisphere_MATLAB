function [x,y,z,d,N_count] = gensph_equidist(r,n_points)

% generate points equidistantly spaced along the surpace of a sphere with
% radius "r", returns [x,y,z] cartesian points
%algorithm based on "How to generate equidistributed points on the surface
%of a sphere" by Markus Deserno

cartcor=[];
N_count=1;
a=(4*pi*r^2)/n_points;
d=sqrt(a);
M_theta=round(pi/d);
d_theta=pi/M_theta;
d_phi=a/d_theta;

for i=(1:M_theta)
    theta = pi*(i-1+0.5)/M_theta; %subtract 1 to adjust for matlab indexing
    M_phi=round(2*pi*sin(theta)/d_phi); 
    for j=(1:M_phi)
        phi=2*pi*(j-1)/M_phi;
        cartcor(N_count,:)=[r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta)];
        N_count=N_count+1;
    end
end

x=cartcor(:,1);
y=cartcor(:,2);
z=cartcor(:,3);

end