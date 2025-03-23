%%%% 3D topology optimization code with skip of DOFs %%%%
% Example: topo3Skip(60,20,10,0.05,0.6,100,1)
function topo3Skip(nelx,nely,nelz,volfrac,rs,nloop,symmdraw)
%% 初始化参数
nelx=60;
nely=20;
nelz=10;
volfrac=0.05;
rs=0.7; 
nloop=100;
symmdraw=1;
penal = 3;                            % Penalization factor 
eps0 = 1e-4;                          % Skip level
eta = 0.6;                            % Damping coefficient, initial value
etamax = 1;                           % and max values
nu = 0.3;                             % Poisson's ratio
%% Load and support fixed dofs
elesize = [nely nelx nelz];
nel = nelx*nely*nelz;
nElemDof = 24;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
% MBB quarter
i = 0; j = nely; k = 0;                                % Coordinates
loadnid = k*(nelx+1)*(nely+1)+i*(nely+1)+(nely+1-j);   % Node IDs
loaddof = 3*loadnid(:) - 1;                            % DOFs
[j,k] = meshgrid(1:nely+1,1:nelz+1);                   % Coordinates
fixednid = (k-1)*(nely+1)*(nelx+1)+j;                  % Node IDs
fixeddof = [3*fixednid(:)-2; 3*fixednid(:)-3];         % DOFs
[i,j] = meshgrid(1:nelx+1,1:nely+1);                   % Coordinates
fixednid = i*(nely+1)+(nely+1-j);                      % Node IDs
fixeddof = [fixeddof(:); 3*fixednid(:)-3];             % DOFs
i = nelx; j = 0; k = 0;
fixednid = k*(nelx+1)*(nely+1)+i*(nely+1)+(nely+1-j);  % Node IDs
fixeddof = [fixeddof(:); 3*fixednid(:)-1];             % DOFs
F = sparse(loaddof,1.,-1.,ndof,1);
%%  有限元分析
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofMat = int32(repmat(3*nodeids(:)+1,1,24)+ ...
	repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
	3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nel,1));

KE = element_stiffness(nu);
U = zeros(ndof,1);
freedof = setdiff(1:ndof,fixeddof);
clearvars loaddof fixeddof nodegrd nodeids;
iK = reshape(kron(edofMat,ones(nElemDof,1,'int32'))',[],1);
jK = reshape(kron(edofMat,ones(1,nElemDof,'int32'))',[],1);
posD = int32(find(jK == iK));             % Positions of diagonal of matrix
iD = iK(posD);
posK = int32(find(jK < iK)); % Positions of lower triangular part of matrix
iK = iK(posK);
jK = jK(posK);
%% 灵敏度分析
function [ce,skip] = sensitivity()
	% Assemble half of global stiffness matrix
	y = x.^penal;	                      % Young modulus
	y(y < eps0*volfrac) = 0;            % is 0 for skipped elements
		
	sK = KE(:)*y'; 
	K = sparse(double(iK),double(jK),sK(posK),ndof,ndof); % Lower triangular
	d = full(sparse(double(iD),1,sK(posD))); % Diagonal
	clearvars sK;
	% Skip DOFs
	skip = find(d == 0);                % Skipped DOFs
  remain = setdiff(freedof,skip);     % Remaining DOFs
	% minimize bandwidth of matrix
	p = remain(symrcm(K(remain,remain)));
	K = K(p,p);
	d = d(p);
	K = K+diag(sparse(d))+K';
	% Preconditioner
	alfa = 0;
	while 1
		try
			P = ichol(K,struct('type','ict','droptol',0.001,'diagcomp',alfa));
			break;
		catch
			alfa = (alfa+0.0001)*2;
		end
	end
	
	% Solve for remaining DOFs only
	[U(p),~] = pcg(K,F(p),[],10000,P,P',U(p));
	ue = U(edofMat);
  ce = max(0,y.*sum(ue*KE.*ue,2));    % Element energy
end
%% Plot faces
function plotFaces()
	clf; caxis([0,1]); colormap hot; set(gcf,'GraphicsSmoothing','off');
	%% v is x with void margins
	v = zeros(elesize+2,'single');
	s1 = elesize+1;
	v(2:s1(1),2:s1(2),2:s1(3)) = reshape(x,elesize);
	%% x0 is cut level to plot only volfrac part
	n = ceil(volfrac*nel);
	s = sort(x,'descend');
  x0 = (s(n) + s(n+1))/2;
	%% Find out number of boundary faces and reserve memory
  s = (v < x0);
	n = nnz(diff(s,1,1)) + nnz(diff(s,1,2)) + nnz(diff(s,1,3));
  vx = zeros(n,4,'int32');     % x coordinates
	vy = zeros(n,4,'int32');     % y coordinates
	vz = zeros(n,4,'int32');     % z coordinates
	col = zeros(n,1,'single');   % colors
	%% function to add faces
  n = 0;
  d = int32([[0 0 0 0]; [0 0 1 1]; [0 1 1 0]; [1 1 1 1];]);
	function addFace(dj,di,dk);
		n = n+1; 
		vx(n,:) = i+d(di,:);
		vy(n,:) = j+d(dj,:);
		vz(n,:) = k+d(dk,:);
		col(n) = v(j,i,k);
  end
  %% Scan neighbors of 'solid' elements to detect boundary faces
  for k = 2:s1(3);  for i = 2:s1(2);  for j = 2:s1(1);
		if ~s(j,i,k)
			if s(j,i,k-1) addFace(2,3,1); end
			if s(j,i,k+1) addFace(2,3,4); end
			if s(j,i-1,k) addFace(2,1,3); end
			if s(j,i+1,k) addFace(2,4,3); end
			if s(j-1,i,k) addFace(1,3,2); end
			if s(j+1,i,k) addFace(4,3,2); end
		end
	end;  end;  end;
	%% Draw
  patch(vx',vz',-vy',col);
	if symmdraw
		patch(vx',4-vz',-vy',col); % symmetric draw for negative y and z
	end;
  light('Position',[1 0 1]); view([30,30]); 
	axis equal; axis tight; axis off; drawnow;
end
x = repmat(volfrac,nel,1);              % Densities


%% 迭代过程
lm = 0;
reta = (etamax/eta)^(1/nloop);          % Augment of damping
for loop = 1:nloop
	%% 高斯密度滤波
	x = imgaussfilt3(reshape(x,elesize),rs); x = x(:);
  %% 灵敏度分析
  [ce,skip] = sensitivity();
  %% Optimality criteria update of densities
	x = x.*(ce.^eta);
	xNew = @(a) min(1,x*exp(a));
	delta = @(a) sum(xNew(a)) - volfrac*nel;
	lm = fzero(delta,lm);
	x = xNew(lm);
	eta = min(etamax,eta*reta);         % New value of damping coefficient
  %% Info
  fprintf('%u\t C: %.6g\t Skip: %u%%\n', ...
    loop,sum(ce),round(100*numel(skip)/numel(freedof)));
	plotFaces();
end
end
%% 计算单元刚度矩阵
function [KE] = element_stiffness(nu)
	a = [6 -3 3 -6 -6 -10 -3 4 3 6 -8 32 -4 -8];
	b = [0 0 0 0 24 12 12 0 -12 -24 12 -48 12 0];
	k = (a + nu*b)/(72*(1+nu)*(1-2*nu));

	i1 = [12 1 1 14 5 5; 1 12 1 10 8 3; 1 1 12 10 3 8;
        14 10 10 12 4 4; 5 8 3 4 12 1; 5 3 8 4 1 12];
	i2 = [6 4 7 8 10 3; 4 6 7 5 14 5; 9 9 13 3 10 8;
        8 5 2 6 1 9; 10 14 5 1 6 7; 2 10 8 7 9 13];
	i3 = [8 3 10 6 7 4; 3 8 10 9 13 9; 5 5 14 4 7 6;
        6 9 1 8 2 5; 7 13 9 2 8 10; 1 7 6 10 5 14];
	i4 = [11 2 2 13 9 9; 2 11 2 7 6 4; 2 2 11 7 4 6;
        13 7 7 11 3 3; 9 6 4 3 11 2; 9 4 6 3 2 11];
	i5 = [12 1 4 14 5 10; 1 12 4 10 8 2; 4 4 12 5 2 8;
        14 10 5 12 4 1; 5 8 2 4 12 4; 10 2 8 1 4 12];
	i6 = [11 2 3 13 9 7; 2 11 3 7 6 1; 3 3 11 9 1 6;
        13 7 9 11 3 2; 9 6 1 3 11 3; 7 1 6 2 3 11];

	KE = k([i1  i2  i3  i4;
					i2' i5  i6  i3';
					i3' i6  i5  i2';
					i4  i3  i2  i1]);
end

