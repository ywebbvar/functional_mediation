function [out] = fMediation_ML(x,y,m,nbasis,norder, varargin)
% function [out] = fMediation_ML(x,y,m)
%
% Performs functional mediation analysis
%
% In the current implementation: x and y are assumed to be vectors and m is timeseries data
%
% Inputs:
%
% x     -   1 X N vector
% y     -   1 X N vector
% m     -   len (time) x N (repetitions) vector of observations 
%
% Outputs: Saved in fields of out.(fieldname)
%
% afun - functional regression coefficient correesponding to a path 
% bfun - functional regression coefficient correesponding to b path 
% ab   - ab effect
% c    - c effect
% cp   - c' effect
%
% By Martin Lindquist, December 2008
%
% EXAMPLES:
% 
% Example 1:
%
% nbasis = 10;
% norder = 6;
% m = zeros(20,60);
% h = HRF(1);
% Z = [zeros(20,1); h'; zeros(10,1)];
% x = unifrnd(0,10,20,1);
% for i=1:20, m(i,:) = (x(i)+normrnd(0,1))*Z' + normrnd(0,1,1,60); end;
% y = sum(m.*repmat(Z',20,1),2) + normrnd(0,1,20,1);
% y = y'; x = x'; m = m';
%
% out = fMediation_ML(x,y,m,nbasis,norder);
% 
%
% Example 2:
% 
% Z = zeros(20,60);
% x = unifrnd(0,10,20,1);
% for i=1:20, 
%     Run = zeros(60,1);
%     Run(21:round(21+x(i))) = 1;
%     s = modifiedconv(1,Run);        
%     Z(i,:) = s'./max(s);
% end;
% 
% m = Z + normrnd(0,0.5,20,60);
% y = sum(m.*Z,2);
% y = y'; x = x'; m = m';
% 
% out = fMediation_ML(x,y,m);


[len, N] = size(m);
T = 1;
timevec = linspace(0,T,len);

% Create bspline basis set
% nbasis = 20;
% norder = 6;
estimate  = 2;
lambda    = 10;  
pen = 0.1;

for varg = 1:length(varargin)
        if ischar(varargin{varg})
            switch varargin{varg}
                % reserved keywords
                case {'penalty'}, pen = varargin{varg+1};
                case {'estimate'}, estimate = varargin{varg+1};
                case {'lambda'} , lambda = varargin{varg+1};
                otherwise, disp(['Unknown input string option: ' varargin{varg}]);
            
            end
        end
end



basis = create_bspline_basis([0,T], nbasis, norder);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  x -> m  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Convert time series data to functional data
mfd = data2fd(m, timevec, basis);
mfdPar = fdPar(mfd);                    % Response variable

% Create Design matrix
conbas = create_constant_basis([0,T]);
confd  = fd(ones(1,N),conbas);

xfdcell = cell(1,2);
xfdcell{1} = confd;             % Constant
xfdcell{2} = x;                 % x

% Create basis set for beta functions
betacell = cell(1,2);
betabasis = create_bspline_basis([0,T],nbasis,norder);
betafd1 = fd(1, conbas);
betacell{1} = fdPar(betafd1);
betafdj = fd(zeros(nbasis,1),betabasis);
betacell{2} = fdPar(betafdj);
% if (estimate > 0 && lambda > 0), 
%     betacell{2} = fdPar(betafdj,estimate,lambda);
% else
%     betacell{2} = fdPar(betafdj);
% end;

% Solve least-squares equation
fRegressCell = fRegress(mfdPar, xfdcell, betacell);
betaestcell = fRegressCell{4};
afun = getfd(betaestcell{2});          % a-function

% % Calculate Standard Error
% errmat = m - eval_fd(timevec,fRegressCell{5});
% Sigma = errmat*errmat'/20;
% DfdPar = fdPar(basis,0,1);
% [fdobj, df, gcv, coef, SSE, penmat, y2cMap] = smooth_basis(timevec,m,DfdPar);
% stderrCell = fRegress_stderr(fRegressCell, y2cMap, Sigma);
% tmp = stderrCell{1};


tfine = linspace(0,T,len)';
af = eval_fd(tfine,afun);               % Evaluate a-function
% a_stderr = eval_fd(timevec,tmp{2});     % Std Error of a-function
a = sum(af)*(tfine(2)-tfine(1));       % Integral of a-function: a = \int af(t) dt
                                       % gives a-path
ResM = (m - eval_fd(tfine,fRegressCell{5}));
                                       
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% x, m -> y %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Create Design matrix
mfdcell = cell(1,3);
mfdcell{1} = confd;
basis = create_bspline_basis([0,T], nbasis, norder);
mfdcell{2} = data2fd(m, timevec, basis);
mfdcell{3} = fd(x,conbas);

% Response variable
yfdPar = y';

% Create basis set for beta functions
betacell = cell(1,3);
betafd1 = fd(1, conbas);
betacell{1} = fdPar(betafd1);
betafdj = fd(zeros(nbasis,1),basis);                          %Penalty term needed here
betafdPar = fdPar(betafdj, estimate, lambda);
betacell{2} = betafdPar;
betacell{3} = fdPar(betafd1);

% Solve least-squares equation
fRegressCell = fRegress(yfdPar, mfdcell, betacell);
betaestcell = fRegressCell{4};
bfun = getfd(betaestcell{2});           % b-function


% Calculate Standard Error
errmat = (y - fRegressCell{5}');
ResY = errmat;
Sigma = errmat*errmat'/20;
DfdPar = fdPar(basis,0,1);
[fdobj, df, gcv, coef, SSE, penmat, y2cMap] = smooth_basis(timevec,m,DfdPar);
stderrCell = fRegress_stderr(fRegressCell, y2cMap, Sigma);
tmp = stderrCell{1};

b_stderr = eval_fd(tfine,tmp{2});     % Std Error of b-function


bf = eval_fd(tfine,bfun);               % Evaluate b-function
b = sum(bf)*(tfine(2)-tfine(1));        % Integral of b-function: b = \int bf(t) dt
                                        % gives b-path

abf = eval_fd(tfine,afun.*bfun);        % ab-functions
ab = sum(abf)*(tfine(2)-tfine(1));      % Integral of ab-function: ab = \int af(t)bf(t) dt
                                        % gives ab-path

                                        
tmp = eval_fd(tfine,getfd(betaestcell{3}));
cp = tmp(1);                            % c'-path


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  x -> y  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Create Design matrix
xfdcell = cell(1,2);
conbas = create_constant_basis([0,T]);
confd  = fd(ones(1,N),conbas);
xfdcell{1} = confd;
xfdcell{2} = fd(x,conbas);

% Response variable
yfdPar = y';

% Create basis set for beta functions
betacell = cell(1,2);
betafd1 = fd(1, conbas);
betacell{1} = fdPar(betafd1);
betacell{2} = fdPar(betafd1);

% Solve least-squares equation
fRegressCell = fRegress(yfdPar, xfdcell, betacell);

betaestcell = fRegressCell{4};


tmp = eval_fd(tfine,getfd(betaestcell{2}));
c = tmp(1);                         % c-path

%IV = c./af;
IV2 = c./af;

A = repmat(af',len,1);
% pen = 0.1;
IV = (inv(A'*A + pen*eye(len) )*A'*ones(len,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Plot Results
% subplot(3,2,1)
% plot(tfine, eval_fd(tfine, afun));
% % hold
% % plot(tfine, af + 2*a_stderr,'g')
% % plot(tfine, af - 2*a_stderr,'g')
% title('a function')
% subplot(3,2,3)
% plot(tfine, eval_fd(tfine,bfun));
% % hold
% % plot(tfine, bf + 2*b_stderr,'g')
% % plot(tfine, bf - 2*b_stderr,'g')
% title('b function')
% subplot(3,2,5)
% plot(tfine,abf);
% title('ab function')
% 
% subplot(2,2,2)
% scatter(x,y)
% title('Scatter plot of x against y'); xlabel('x'); ylabel('y');
% subplot(2,2,4)
% plot(m)
% title('m');
% 
% % Print Results
% fprintf('\n');
% fprintf('Parameter \n');
% fprintf('c        %3.3f \n', c)
% fprintf('cp       %3.3f \n', cp)
% fprintf('ab       %3.3f \n', ab)
% fprintf('\n');
% fprintf('cp + ab  %3.3f \n', cp+ab)

% Output
out = struct('afunction', af, 'a', a, 'bfunction', bf, 'abfunction', abf, 'b', b, 'ab', ab, 'c', c, 'cp', c,'x',x, 'y', y, 'm', m,'tfine',tfine,'b_stderr',b_stderr,'IV',IV,'IV2',IV2,'ResM',ResM,'ResY', ResY);

end