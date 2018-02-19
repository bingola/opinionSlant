function [ahat]=inference(user,graph)
A=graph.Adj;
Nuser=find(A(:,user));
delmuf=0;
delAuvf=0;
for sample=1:1
  USS=graph.USSx{sample};
  T=max(USS{user});
  for j=1:length(Nuser)
      if ~isempty(USS{Nuser(j)})
       T=max(T,max(USS{Nuser(j)}));
     end
 end

 g=graph.g{user};
 G=graph.G{user};
 if isempty(g)
  ahat=zeros(length(Nuser)+1,1);
 
 else
      
 L=10;
 %%%set error%%%
 e=5e-2;
 iter=1;
  a(:,1)=ones(length(Nuser)+1,1);
 
  I=(eye(length(Nuser)+1));
   Hk=I;
  alphamin=0.0001;
   alphabb=0.5*rand+0.0001;
  dk=10;
  nux=1e-3;
  Bk=I;
 
 while norm(dk)>e
%       L=norm(dk)
 mu1=a(1,iter);
 a1=a(2:end,iter);

 mx=mu1+g'*a1;
 yvec=repmat(mx',length(Nuser),1);
 ymat=g./yvec;
 delAuv=G'-(sum(ymat'))';
 delmu=T-sum(1./(mu1+g'*a1));
 wk=[delmu;delAuv];
 
 alphak=min(alphamin,max(alphamin,alphabb));
 dk=projected(a(:,iter)-alphak*wk)-a(:,iter);
 Dqk=wk;
 alpha=1;
%  Bk=inv(Hk);
 while(alpha*dk'*Dqk+alpha^2*dk'*Bk*dk>nux*alpha*Dqk'*dk&& alpha>0)
  alpha=alpha-0.1;
 end
 
 a(:,iter+1)=a(:,iter)+alpha*dk;
 sk=a(:,iter+1)-a(:,iter);
 iter=iter+1;
 mu1=a(1,iter);
 a1=a(2:end,iter);
 mx=mu1+g'*a1;
 yvec=repmat(mx',length(Nuser),1);
 
 
  yk=Bk*sk;
  Hk=(I-sk*yk'/(yk'*sk))*Hk*((I-sk*yk'/(yk'*sk)))'+sk*sk'/(yk'*sk);

 
 alphabb=yk'*yk/(sk'*yk);
  
 Bk=Bk-Bk*(sk*sk')*Bk/(sk'*Bk*sk)+yk*yk'/(yk'*sk);
 end

 ahat=a(:,iter);
 end
end
end