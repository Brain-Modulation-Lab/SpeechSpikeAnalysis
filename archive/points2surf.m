function s = points2surf( Points, dx, dy )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

x=Points(:,1);
y=Points(:,2);
z=Points(:,3);


x_edge=[floor(min(x)):dx:ceil(max(x))];
y_edge=[floor(min(y)):dy:ceil(max(y))];
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x,y,z,X,Y);
% The following line of code is if you use JE's gridfit:
% Z=gridfit(x,y,z,x_edge,y_edge);

%NOW surf and mesh will work...

s = mesh(X,Y,Z);
%mesh(X,Y,Z)
end

