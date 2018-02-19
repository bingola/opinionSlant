function graphout=simulateHawkesEventsOnly(graphx,graphType,averageEvents)
clc
dR=rand(1000,1);
 for sampleNo=1:1  
N=length(graphx.mu);
b=graphx.B;
mu=graphx.mu;
avg=averageEvents;
MaxUserEvents=avg*N;
t0=0; %initial instant

% Initialization
USS=cell(N,1);
tu=zeros(1,N);%number of events per user
lambdaTransient=zeros(1,N);
lambdaAtTimeEqualS=zeros(1,N);

for n=1:N
    USS{n}=zeros(MaxUserEvents,1);    
end
numUses=MaxUserEvents;
tu_events=zeros(2,numUses); % tuple with the events (time,user)


%% first event generation

lambdaIni=mu;

t=1;
u=rand(1);
s=-log(u)/sum(lambdaIni);
d=dR(sampleNo);


I=sum(lambdaIni);
for i=1:N
    if i==1 && d<= sum(lambdaIni(1:i))/I
       user=1;
    elseif  i>1 && d>(sum(lambdaIni(1:i-1))/I) && d<=sum(lambdaIni(1:i))/I
        user=ceil(i);
    end
end

USS{user}(1,1)=s;
tu(user)=tu(user)+1;
tu_events(1,1)=s; %% time
tu_events(2,1)=user; %% user
nuVal=2;
%% general subroutine
 while sum(tu)<MaxUserEvents
%      sum(tu)
% while min(tu-minNumberPerUser)<0
     % updating intensities with the last event
    lambdaTransient=lambdaTransient+b(user,:);
    
    
    %sampling a new event
    lambda2=lambdaIni+lambdaTransient;
    lambda2((lambdaTransient<0))=0;
    
    I=sum(lambda2);lambdaAtTimeEqualS
    u=rand(1);
    s=s-log(u)/I;

    lambdaAtTimeEqualS=lambdaIni+lambdaTransient*exp(nuVal*log(u)/I);
    Is=sum(lambdaAtTimeEqualS);
    lambdaTransient=lambdaTransient*exp(nuVal*log(u)/I);

    d=rand(1);
    if d<= Is/I
        for i=1:N
            if i==1 && d<= sum(lambdaAtTimeEqualS(1:i))/I
              user=1;
            elseif  i>1 && d>(sum(lambdaAtTimeEqualS(1:i-1))/I) && d<=sum(lambdaAtTimeEqualS(1:i))/I
              user=i;
            end
        end
        tu(user)=tu(user)+1;
        USS{user}(tu(user),1)=s;
    else
     
     while d> Is/I
        lambda2=lambdaIni+lambdaTransient;
        lambda2((lambdaTransient<0))=0;
        I=sum(lambda2);
        u=rand(1);
        s=s-log(u)/I;
        lambdaAtTimeEqualS=lambdaIni+lambdaTransient*exp(nuVal*log(u)/I);
        Is=sum(lambdaAtTimeEqualS);
        lambdaTransient=lambdaTransient*exp(nuVal*log(u)/I);
        d=rand(1);
       
        for i=1:N
            if i==1 && d<= sum(lambdaAtTimeEqualS(1:i))/I
              user=1;
            elseif  i>1 && d>(sum(lambdaAtTimeEqualS(1:i-1))/I) && d<=sum(lambdaAtTimeEqualS(1:i))/I
              user=i;
            end
        end
%        lambdaTransient=lambdaTransient*exp(nuVal*log(u)/I);
     end
        tu(user)=tu(user)+1;
        USS{user}(tu(user),1)=s;
    end
        tu_events(1,t)=s;
        tu_events(2,t)=user;  
end


% events per user
for n=1:N
    USS{n}=USS{n}(1:tu(n),:);
end
if sampleNo==1
graphout=graphx;
end
graphout.USSx{sampleNo}=USS;
graphout.tux{sampleNo}=tu;


end
%%
clear g
clear G
g{N}=[];
G{N}=[];
clc
w=nuVal;
A=graphout.Adj;
for user=1:N
  if ~isempty(USS{user})  
   Nuser=find(A(:,user));
   T=max(USS{user});
   for j=1:length(Nuser)
      if ~isempty(USS{Nuser(j)})
       T=max(T,max(USS{Nuser(j)}));
     end
   end
  
 
   for j=1:length(Nuser) 
      
    g1=(ones(length(USS{Nuser(j)}),1)*USS{user}'-USS{Nuser(j)}*ones(1,length(USS{user})));
    g1(g1<0)=Inf;
%     if user==151
%       size(g1)
%     end
    g{user}(j,:)=sum(exp(-g1),1);
    tj=USS{Nuser(j)};
    K=sum((1-exp(-w*(T-tj))/w));   
    G{user}(j)=K;
   end
  end
end





graphout.g=g;
graphout.G=G;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=(1e-3)*T;
ts=0:dt:T;
clear g
clear G
g{N}=[];
G{N}=[];
clc
w=nuVal;
A=graphout.Adj;
for user=1:N
%   if ~isempty(USS{user})  
%    Nuser=find(A(:,user));
%    T=max(USS{user});
%    for j=1:length(Nuser)
%       if ~isempty(USS{Nuser(j)})
%        T=max(T,max(USS{Nuser(j)}));
%      end
%    end
  
 
   for j=1:length(Nuser) 
      
    g1=(ones(length(USS{Nuser(j)}),1)*ts-USS{Nuser(j)}*ones(1,length(ts)));
    g1(g1<0)=Inf;
%     if user==151
%       size(g1)
%     end

     g{user}(j,:)=sum(exp(-g1),1)';
%     tj=USS{Nuser(j)};
%     K=sum((1-exp(-w*(T-tj))/w));   
%     G{user}(j)=K;
   end
%    size(b(Nuser,user)'*g{user})
   
   lambdastar(user,:)=graphx.mu(user)*ones(length(ts),1)'+b(Nuser,user)'*g{user};
   
  
  
end
graphout.lambdastar=lambdastar;


%%%%%%%%%%%%%%%%%%%%%%%



FILE=strcat('../data/',graphType,'512Graphout');
avgm=avg;
if avg<1
    avgm=ceil(1/avg)+1;
end
eval(['save -v7.3 ',FILE,num2str(avgm),' graphout;']);

