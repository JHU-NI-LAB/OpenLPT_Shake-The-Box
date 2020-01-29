function d = point_to_line(pt, v1, v2)
%       a = v1 - v2;
%       b = pt - repmat(v2,[size(pt,1) 1]);
%       aa=repmat(a, [size(pt,1) 1]);
%       aa=[aa,zeros(size(pt,1),1)];
%       bb=[b,zeros(size(pt,1),1)];
%       d = norm(cross(aa,bb)) / norm(aa);

for i=1:size(pt,1)
    a=pt(i,:);
    d(i,1)=abs( det([a-v1;v2-v1]) )/norm(v2-v1);
end

