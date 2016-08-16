addpath ../fem_util/
addpath ../C_files/
% addpath ../data/
addpath ../meshing/
% addpath ../post-processing/
% addpath ../fem-functions/
% addpath ../analytical-solutions/
addpath ../nurbs-util/
close all
clc
clear all
close all
tic
global p
p       = 1
nctrPts=160;
% nxi = 19;
nxi = nctrPts+1-p;
knotVec=[zeros(1,p) linspace(0,1,nxi) ones(1,p)];
controlPts = [linspace(0,1,nxi+p-1);zeros(1,nxi+p-1)]';
nctrPts=size(controlPts)
% p       = 4;
% % knotVec    = [0 0 0 0 1 1 1 1];
% knotVec = [zeros(1,p+1) ones(1,p+1)]
% controlPts = [linspace(0,1,p+1);zeros(1,p+1)]'
% controlPts = [0 0;
%               1/3 0; 
%               2/3 0;
% %               1/2 0;
%               1 0];
% controlPts = [0 0;
%               1/4 0; 
%               1/2 0;
%               3/4 0;
%               1 0];
          
noGPs   = p+1; % number of Gauss points

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));
refinement =0;
num_max_iter = 1; % Max. number of summands for the approach
TOL = 1.0E-8; 
%%%% INITIAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[];Y=[];
generateIGA1DMesh;
%===========plot control mesh===============
noPts      = 2;
xi         = linspace(0,max(knotVec),noPts);
sCurve     = zeros(2,noPts);

for i=1:noPts
  sCurve(1,i) = NURBSinterpolation(xi(i), p, knotVec, controlPts(:,1), controlPts(:,3));
  sCurve(2,i) = NURBSinterpolation(xi(i), p, knotVec, controlPts(:,2), controlPts(:,3));
end

xi = unique(knotVec);
for i=1:numel(xi)
  knot(1,i) = NURBSinterpolation(xi(i), p, knotVec, controlPts(:,1), controlPts(:,3));
  knot(2,i) = NURBSinterpolation(xi(i), p, knotVec, controlPts(:,2), controlPts(:,3));
end
figure;
plot(sCurve(1,:),sCurve(2,:),'b-','LineWidth',1.8);
plot(knot(1,:),knot(2,:),'s');
plot(controlPts(:,1),controlPts(:,2),'r-o',...
                'MarkerEdgeColor','y',...
                'MarkerFaceColor','b',...
                'MarkerSize',10);
noCtrPts = size(controlPts,1);% no of control points
noElems  = size(elConn,1);    % no of elements

% initialization
M = zeros(noCtrPts,noCtrPts); % global stiffness matrix 
dM = zeros(noCtrPts,noCtrPts); % global stiffness matrix 
% u = zeros(noCtrPts,1);        % displacement vector
% f = zeros(noCtrPts,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(noGPs, 'GAUSS', 1 ); % 2 point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

% disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)
for e=1:noElems
   xiE   = elRange(e,:); % [xi_i,xi_i+1]
   conn  = elConn(e,:);
   noFns = length(conn);
 
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));% coord in parameter space  
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      [N,dNdxi] = NURBS1DBasisDers(Xi,p,knotVec,controlPts(:,3)');

      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      jacob1 = dNdxi*controlPts(conn,1:2);
      J1     = norm (jacob1);
      dNdx   = (1/J1)*dNdxi;
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      M(conn,conn) = M(conn,conn) + N' * N* J1 * J2 * wt;
      
      % compute the external force, kind of body force
%       X       = N * controlPts(conn,1:2);
%       bx      = X(1); % f(x)=x
      dM(conn,conn) = dM(conn,conn) + dNdx' * dNdx * J1 * J2 * wt; 
    end
end
% Ah=ones(1,noCtrPts);
% Bh=ones(1,noCtrPts);
x=controlPts(:,1)';
Ah=sin(pi*x);
Bh=sin(pi*x);
Ax=M*Ah'; 
By=M*Bh';
%%% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC_R=(2:1:noCtrPts-1)';
CC_S=(2:1:noCtrPts-1)';
%%%% PGD SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_iter = 0; iter = zeros(1); Aprt = 0; Error_iter = 1.0; 
while Error_iter>TOL && num_iter<num_max_iter
    num_iter = num_iter + 1;
    R0=ones(noCtrPts,1);R0(1,1)=0;R0(end,1)=0;
    S0=ones(noCtrPts,1);S0(1,1)=0;S0(end,1)=0;
%%% ENRICHMENT STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R,S,iter(num_iter)] = enrich_step(M,dM,M,dM,Ax,By,num_iter,X,Y,R0,S0,CC_R,CC_S,TOL,controlPts(:,1),controlPts(:,1));
    X(:,num_iter)=R;
    Y(:,num_iter)=S;
%%% ERROR STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Error_iter=norm(X(:,num_iter))*norm(Y(:,num_iter));
    Aprt = max(Aprt,sqrt(Error_iter));
    Error_iter = sqrt(Error_iter)/Aprt;
%     fprintf(1,'%dst summand in %d iterations with a weight of %f\n',...
%         num_iter,iter(num_iter),sqrt(Error_iter));
end
u=Y*X';

% xx=linspace(0,1,99);
xx=controlPts(:,1)';
X_ex=sqrt(1/(2*pi^2))*(sin(pi*xx))';
Y_ex=sqrt(1/(2*pi^2))*(sin(pi*xx))';
u_ex=1/(2*pi^2)*(sin(pi*xx))'*sin(pi*xx);
e=norm(u_ex-u)/norm(u_ex);
e1=sqrt((X'*M*X)*(Y'*M*Y) + (X_ex'*M*X_ex)*(Y_ex'*M*Y_ex) - 2*(X'*M*X_ex)*(Y'*M*Y_ex));
e2=sqrt((X'*M*X)*(Y'*M*Y) + (X_ex'*M*X_ex)*(Y_ex'*M*Y_ex) - 2*(X'*M*X_ex)*(Y'*M*Y_ex))/sqrt((X_ex'*M*X_ex)*(Y_ex'*M*Y_ex))
% e2=norm(u_ex-u)
% save('p.mat','u','X','Y','x');
figure()
surf(x,x,Y*X');
toc
figure()
surf(xx,xx,u_ex)
% saiso=norm(u_ex-Y*X')/norm(u_ex)

