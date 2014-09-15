tic

iters = 500;
h = HRF(1);
Z = [zeros(20,1); h'; zeros(10,1)];
sub = 20;
nbasis = 30;
norder = 6;             % Set order of spline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%P = zeros(iters,1);
Pall = zeros(iters,length(Z));
Pall2 = zeros(iters,length(Z));
Flag = zeros(iters,1);

C = zeros(iters,1);
CP = zeros(iters,1);
A = zeros(iters,length(Z));
B = zeros(iters,length(Z));


for i=1:iters

    % Null Data
    [a,b]=sort(rand(20,1));
    X = zeros(20,1);
    X(b(1:10)) = 2;
    X = X - 1;
    M = normrnd(0,1,sub,length(Z));
    Y = X + normrnd(0,1,sub,1);

    try
        out = fMediation_ML(X',Y',M',nbasis,norder);       % Perform functional mediation
        stats = fbootstrap_ML(out, 1000,nbasis,norder);         % Perform stats
        
        Pall(i,:) = 2*min(mean(stats.abfMat < 2*repmat(out.abfunction',1000,1)), mean(stats.abfMat > 2*repmat(out.abfunction',1000,1)));        
        Pall2(i,:) = 2*min(mean(stats.abfMat < repmat(out.abfunction',1000,1)), mean(stats.abfMat > repmat(out.abfunction',1000,1)));        
        A(i,:) = out.afunction';
        B(i,:) = out.bfunction';
        C(i) = out.c;
        CP(i) = out.cp;
    catch
        Flag(i) = 1;
    end
    
    i
    
end

save 'Sim_1_rev_R2' A B C CP Pall Pall2 Flag


%%

Pall = zeros(iters,length(Z));
Flag = zeros(iters,1);

C = zeros(iters,1);
CP = zeros(iters,1);
A = zeros(iters,length(Z));
B = zeros(iters,length(Z));


for i=1:iters

    % Null Data
    [a,b]=sort(rand(20,1));
    X = zeros(20,1);
    X(b(1:10)) = 2;
    X = X - 1;
%   alp = normrnd(0,0.1);
%   M = (X+1+alp)*Z' + normrnd(0,1,sub,length(Z));
    M = (X+1)*Z' + normrnd(0,1,sub,length(Z));
    Y = X + normrnd(0,1,sub,1);

    try
        out = fMediation_ML(X',Y',M',nbasis,norder);       % Perform functional mediation
        stats = fbootstrap_ML(out, 1000,nbasis,norder);         % Perform stats

        Pall(i,:) = 2*min(mean(stats.abfMat < 2*repmat(out.abfunction',1000,1)), mean(stats.abfMat > 2*repmat(out.abfunction',1000,1)));        
        Pall2(i,:) = 2*min(mean(stats.abfMat < repmat(out.abfunction',1000,1)), mean(stats.abfMat > repmat(out.abfunction',1000,1)));        
        A(i,:) = out.afunction';
        B(i,:) = out.bfunction';
        C(i) = out.c;
        CP(i) = out.cp;
        
    catch
        Flag(i) = 1;
    end
    
    i
    
end

save 'Sim_2_rev_R2' A B C CP Pall Pall2 Flag


%%

Pall = zeros(iters,length(Z));
Flag = zeros(iters,1);

C = zeros(iters,1);
CP = zeros(iters,1);
A = zeros(iters,length(Z));
B = zeros(iters,length(Z));

for i=1:iters

    [a,b]=sort(rand(20,1));
    X = zeros(20,1);
    X(b(1:10)) = 2;
    X = X - 1;
    M = (zeros(20,1)+1)*Z' + normrnd(0,1,sub,length(Z));
    Y = sum(M(:,21:40),2) + normrnd(0,1,sub,1);

    try
        out = fMediation_ML(X',Y',M',nbasis,norder);       % Perform functional mediation
        stats = fbootstrap_ML(out, 1000,nbasis,norder);         % Perform stats

        Pall(i,:) = 2*min(mean(stats.abfMat < 2*repmat(out.abfunction',1000,1)), mean(stats.abfMat > 2*repmat(out.abfunction',1000,1)));        
        Pall2(i,:) = 2*min(mean(stats.abfMat < repmat(out.abfunction',1000,1)), mean(stats.abfMat > repmat(out.abfunction',1000,1)));        
        A(i,:) = out.afunction';
        B(i,:) = out.bfunction';
        C(i) = out.c;
        CP(i) = out.cp;

    catch
        Flag(i) = 1;
    end

    i
end

save 'Sim_3_rev_R2' A B C CP Pall Pall2 Flag


%%
Pall = zeros(iters,length(Z));
Flag = zeros(iters,1);

C = zeros(iters,1);
CP = zeros(iters,1);
A = zeros(iters,length(Z));
B = zeros(iters,length(Z));

for i=1:iters

    [a,b]=sort(rand(20,1));
    X = zeros(20,1);
    X(b(1:10)) = 2;
    X = X - 1;    
    M = (X+1)*Z' + normrnd(0,1,sub,length(Z));
    Y = sum(M(:,21:40),2) + normrnd(0,1,sub,1);

    try
        out = fMediation_ML(X',Y',M',nbasis,norder);       % Perform functional mediation
        stats = fbootstrap_ML(out, 1000,nbasis,norder);         % Perform stats
        
        Pall(i,:) = 2*min(mean(stats.abfMat < 2*repmat(out.abfunction',1000,1)), mean(stats.abfMat > 2*repmat(out.abfunction',1000,1)));        
        Pall2(i,:) = 2*min(mean(stats.abfMat < repmat(out.abfunction',1000,1)), mean(stats.abfMat > repmat(out.abfunction',1000,1)));        
        A(i,:) = out.afunction';
        B(i,:) = out.bfunction';
        C(i) = out.c;
        CP(i) = out.cp;

    catch
        Flag(i) = 1;
    end

    i
end

save 'Sim_4_rev_R2' A B C CP Pall Pall2 Flag


%%
iters = 500;
P = zeros(iters,1);
Pall = zeros(iters,length(Z));
Pall_IV = zeros(iters,length(Z));
Flag = zeros(iters,1);
IV = zeros(iters,length(Z));
abIV = zeros(iters,length(Z));

nbasis = 50;


for i=1:iters

    [a,b]=sort(rand(20,1)); 
    X = zeros(20,1);
    X(b(1:10)) = 2;
    X = X - 1;

    eps = normrnd(0,1,sub,length(Z));
    M = X*Z' + eps;
    eta = sum(eps(:,21:30),2) + normrnd(0,1,sub,1);
    Y = sum(M(:,21:30),2) + eta;

    try
      out = fMediation_ML(X',Y',M',nbasis,norder);            % Perform functional mediation
      IV(i,:) = out.IV;
      ab_IV = out.IV.*out.afunction;
      stats = fbootstrap_ML(out, 1000,nbasis,norder);         % Perform stats
      abIV(i,:) = 2*min(mean(stats.abIVMat < 2*repmat(out.ab_IV',1000,1)), mean(stats.abIVMat > 2*repmat(out.ab_IV',1000,1)));   
      Pall_IV(i,:) = 2*min(mean(stats.IVMat < 2*repmat(out.IV',1000,1)), mean(stats.IVMat > 2*repmat(out.IV',1000,1)));  
      Pall(i,:) = 2*min(mean(stats.abfMat < 2*repmat(out.abfunction',1000,1)), mean(stats.abfMat > 2*repmat(out.abfunction',1000,1)));   
      Pall2(i,:) = 2*min(mean(stats.abfMat < repmat(out.abfunction',1000,1)), mean(stats.abfMat > repmat(out.abfunction',1000,1)));        

    catch
        Flag(i) = 1;
    end

    i
end

save 'Sim_6_rev_R2' P Pall Pall2 IV Pall_IV Flag

toc
