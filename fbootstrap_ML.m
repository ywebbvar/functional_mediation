function [stats] = fbootstrap_ML(out, reps,nbasis,norder, varargin)

m = out.m;
x = out.x;
y = out.y;

[len, N] = size(m);
T = 1;
timevec = linspace(0,T,len);
    
% Create bspline basis set
% nbasis = 20;
% norder = 6;

estimate  = 2;
lambda    = 10;  

for varg = 1:length(varargin)
        if ischar(varargin{varg})
            switch varargin{varg}
                % reserved keywords
                case {'estimate'}, estimate = varargin{varg+1};
                case {'lambda'} , lambda = varargin{varg+1};
                otherwise, disp(['Unknown input string option: ' varargin{varg}]);
            
            end
        end
end

basis = create_bspline_basis([0,T], nbasis, norder);    

abMat = zeros(reps,1);
abfMat = zeros(reps,len);
afMat = zeros(reps,len);
IVMat = zeros(reps,len);
abIVMat = zeros(reps,len);

for i=1:reps,

    samp = ceil(unifrnd(0,N,N,1));
    mm = m(:,samp);
    xx = x(samp);
    yy = y(samp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%  x -> m  %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Convert time series data to functional data
    mfd = data2fd(mm, timevec, basis);
    mfdPar = fdPar(mfd);                    % Response variable

    % Create Design matrix
    conbas = create_constant_basis([0,T]);
    confd  = fd(ones(1,N),conbas);
    xfdcell = cell(1,2);
    xfdcell{1} = confd;             % Constant
    xfdcell{2} = xx;                 % x

    % Create basis set for beta functions
    betacell = cell(1,2);
    betabasis = create_bspline_basis([0,T],nbasis,norder);
    betafd1 = fd(1, conbas);
    betacell{1} = fdPar(betafd1);
    betafdj = fd(zeros(nbasis,1),betabasis);
    betacell{2} = fdPar(betafdj);
%    betacell{2} = fdPar(betafdj,estimate,lambda);
    
    % Solve least-squares equation
    fRegressCell = fRegress(mfdPar, xfdcell, betacell);
    betaestcell = fRegressCell{4};
    afun = getfd(betaestcell{2});          % a-function
    tfine = linspace(0,T,len)';
    af = eval_fd(tfine,afun);               % Evaluate a-function
    a = sum(af)*(tfine(2)-tfine(1));       % Integral of a-function: a = \int af(t) dt
                                           % gives a-path

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% x, m -> y %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create Design matrix
    mfdcell = cell(1,3);
    mfdcell{1} = confd;
    basis = create_bspline_basis([0,T], nbasis, norder);
    mfdcell{2} = data2fd(mm, timevec, basis);
    mfdcell{3} = fd(xx,conbas);

    % Response variable
    yfdPar = yy';

    % Create basis set for beta functions
    betacell = cell(1,3);
    betafd1 = fd(1, conbas);
    betacell{1} = fdPar(betafd1);
    betafdj = fd(zeros(nbasis,1),basis);
    betafdPar = fdPar(betafdj, estimate, lambda);
    betacell{2} = betafdPar;
    betacell{3} = fdPar(betafd1);

    % Solve least-squares equation
    fRegressCell = fRegress(yfdPar, mfdcell, betacell);
    betaestcell = fRegressCell{4};
    bfun = getfd(betaestcell{2});           % b-function


    bf = eval_fd(tfine,bfun);               % Evaluate b-function
    b = sum(bf)*(tfine(2)-tfine(1));        % Integral of b-function: b = \int bf(t) dt
                                            % gives b-path




    abf = eval_fd(tfine,afun.*bfun);        % ab-functions
    ab = sum(abf)*(tfine(2)-tfine(1));      % Integral of ab-function: ab = \int af(t)bf(t) dt
                                            % gives ab-path

                                            
    abMat(i) = ab;
    abfMat(i,:) = abf';
    
    
    X = zeros(N,2);
    X(:,1) = 1;
    X(:,2) = x';
    b = regress(y',X);
    IVMat2(i,:) = b(2)./af;
    
    A = repmat(af',len,1);
    pen = 0.1;
    IVMat(i,:) = (inv(A'*A + pen*eye(len) )*A'*ones(len,1));
    abIVMat(i,:) = af'.*IVMat(i,:);

    afMat(i,:) = af';
end                      
    
alpha = 0.05;
down = round(reps*alpha/2);
up = reps-down;
p = mean(abMat < 0);
s = sort(abfMat,1);
abf_lower = s(down,:);
abf_upper = s(up,:);
               
stats = struct('abfunction', out.abfunction, 'abf_upper', abf_upper, 'abf_lower', abf_lower,'pvalue', p,...
    'abMat', abMat, 'abfMat', abfMat, 'afMat', afMat, 'abIVMat', abIVMat, 'reps', reps,'IVMat',IVMat,'IVMat2',IVMat2);


end                                        