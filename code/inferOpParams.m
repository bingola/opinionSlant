function [a]=inferOpParams(user,graph,regVal)
g=graph.gOpinion;
[r,c]=size(graph.OSS{user});
g1=[ones(1,r);g{user}];
 
g=g1;
Y=graph.OSS{user};
%%%Y=Xb%%%%%
X=g';
% [Q,R]=qr(X);
R=X;
L=R'*R;
I=eye(max(size(L)));
a=(regVal*I+L)\(X'*Y);
 

end

