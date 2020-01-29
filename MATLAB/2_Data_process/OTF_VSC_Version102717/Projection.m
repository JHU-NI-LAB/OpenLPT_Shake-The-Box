function z = Projection(par,xy)
% compute 2D gaussian
xx=(xy(1)-par(5)).*cos(par(2))+(xy(2)-par(6)).*sin(par(2));
yy=-(xy(1)-par(5)).*sin(par(2))+(xy(2)-par(6)).*cos(par(2));
z = round(par(1)*exp(-(xx.^2 * par(3) + yy.^2 * par(4))));

end

