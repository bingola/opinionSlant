function a=projected(v)
% a=v;
v(v<0)=0;
v(v>1)=1;
a=v;
end