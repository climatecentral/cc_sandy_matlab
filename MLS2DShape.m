function [PHI, DPHIx, DPHIy] = MLS2DShape(m, nnodes, xI,yI, npoints, x,y, dmI, type, para)
% SHAPE FUNCTION OF 2D MLS APPROXIMATION
%
% SYNTAX: [PHI, DPHI, DDPHI] = MLS1DShape(m, nnodes, xI,yI, npoints, xi,yi, dmI, type, para)
%
% INPUT PARAMETERS
%    m - Total number of basis functions (1: Constant basis;  2: Linear basis;  3: Quadratic basis)
%    nnodes  - Total number of nodes used to construct MLS approximation
%    npoints - Total number of points whose MLS shape function to be evaluated
%    xI,yI(nnodes) - Coordinates of nodes used to construct MLS approximation
%    xi,yi(npoints) - Coordinates of points whose MLS shape function to be evaluated
%    dm(nnodes) - Radius of support of nodes
%    wtype - Type of weight function
%    para  - Weight function parameter
%
% OUTPUT PARAMETERS
%    PHI   - MLS Shpae function
%    DPHIx  - First order derivatives of MLS Shpae function to x
%    DPHIy - First order derivatives of MLS Shpae function to y
%
% INITIALIZE WEIGHT FUNCTION MATRICES
DmI=[];


% INITIALIZE SHAPE FUNCTION MATRICES
PHI   = nan * zeros(npoints, nnodes);
DPHIx  = zeros(npoints, nnodes);
DPHIy = zeros(npoints, nnodes);

parfor_progress(npoints);
% LOOP OVER ALL EVALUATION POINTS TO CALCULATE VALUE OF SHAPE FUNCTION Fi(X)
parfor j = 1 : npoints'
    %j
    
        
    parfor_progress;
    
    
    if isnan(x(j)) || isnan(y(j))
       continue
    end
    
    
    DmI = dmI;
    
    wI   = zeros (1, nnodes);  % Weight funciton
    dwdxI  = zeros (1, nnodes);
    dwdyI = zeros (1, nnodes);
    xII = zeros(1,nnodes);
    yII = zeros(1,nnodes);

	% DETERMINE WEIGHT FUNCTIONS AND THEIR DERIVATIVES AT EVERY NODE
	for i = 1 : nnodes
		[wI(i), dwdxI(i), dwdyI(i)] =rectangleWeight(type, para, x(j),y(j),xI(i),yI(i),DmI(i),DmI(i));
       xII(1,i)=xI(i);
       yII(1,i)=yI(i);
    end
   
   % EVALUATE BASIS p, B MATRIX AND THEIR DERIVATIVES
   if (m == 1)  % Shepard function
      p = [ones(1, nnodes)]; 
      pxy   = [1];
      dpdx  = [0];
      dpdy = [0];
      
      B    = p .* [wI];
      DBdx   = p .* [dwdxI];
      DBdy  = p .* [dwdyI];
   elseif (m == 3)
      p = [ones(1, nnodes); xII;yII]; 
      pxy   = [1; x(j);y(j)];
      dpdx  = [0; 1;0];
      dpdy = [0; 0;1];
      
      B    = p .* [wI; wI;wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI];
   elseif (m == 6)
      p = [ones(1, nnodes); xII;yII; xII.*xII;xII.*yII;yII.*yII]; 
      pxy   = [1; x(j); y(j);x(j)*x(j);x(j)*y(j);y(j)*y(j)];
      dpdx  = [0; 1; 0;2*x(j);y(j);0];
      dpdy = [0;0;1;0;x(j);2*y(j)];
      
      B    = p .* [wI; wI; wI; wI; wI; wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI;dwdxI; dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI;dwdyI; dwdyI;dwdyI];
   else
      error('Invalid order of basis.');
   end
   
   % EVALUATE MATRICES A AND ITS DERIVATIVES
	A   = zeros (m, m);
	DAdx  = zeros (m, m);
	DAdy= zeros (m, m);
	for i = 1 : nnodes
      pp = p(:,i) * p(:,i)';
      A   = A   + wI(i) * pp;
      DAdx  = DAdx  + dwdxI(i) * pp;
      DAdy = DAdy + dwdyI(i) * pp;
    end
   ARcond = rcond(A);
   
   
   while ARcond<=9.999999e-015  %�ж�������
       DmI=1.1*DmI;
 for i = 1 : nnodes
		[wI(i), dwdxI(i), dwdyI(i)] =rectangleWeight(type, para, x(j),y(j),xI(i),yI(i),DmI(i),DmI(i));
       xII(1,i)=xI(i);
       yII(1,i)=yI(i);
 end
   
   % EVALUATE BASIS p, B MATRIX AND THEIR DERIVATIVES
   if (m == 1)  % Shepard function
      p = [ones(1, nnodes)]; 
      pxy   = [1];
      dpdx  = [0];
      dpdy = [0];
      
      B    = p .* [wI];
      DBdx   = p .* [dwdxI];
      DBdy  = p .* [dwdyI];
   elseif (m == 3)
      p = [ones(1, nnodes); xII;yII]; 
      pxy   = [1; x(j);y(j)];
      dpdx  = [0; 1;0];
      dpdy = [0; 0;1];
      
      B    = p .* [wI; wI;wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI];
   elseif (m == 6)
      p = [ones(1, nnodes); xII;yII; xII.*xII;xII.*yII;yII.*yII]; 
      pxy   = [1; x(j); y(j);x(j)*x(j);x(j)*y(j);y(j)*y(j)];
      dpdx  = [0; 1; 0;2*x(j);y(j);0];
      dpdy = [0;0;1;0;x(j);2*y(j)];
      
      B    = p .* [wI; wI; wI; wI; wI; wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI;dwdxI; dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI;dwdyI; dwdyI;dwdyI];
   else
      error('Invalid order of basis.');
   end
   
   % EVALUATE MATRICES A AND ITS DERIVATIVES
	A   = zeros (m, m);
	DAdx  = zeros (m, m);
	DAdy= zeros (m, m);
	for i = 1 : nnodes
      pp = p(:,i) * p(:,i)';
      A   = A   + wI(i) * pp;
      DAdx  = DAdx  + dwdxI(i) * pp;
      DAdy = DAdy + dwdyI(i) * pp;
    end
   ARcond = rcond(A);   
   end   %�ж�������   
   
   
   
   AInv = inv(A);
      
   rxy  = AInv * pxy;
   PHI(j,:) = rxy' * B;   % shape function
    
   drdx  = AInv * (dpdx -DAdx* rxy);
   DPHIx(j,:) = drdx' * B + rxy' * DBdx;   % first order derivatives of shape function with respect to x
   
     drdy = AInv * (dpdy -DAdy* rxy);
   DPHIy(j,:) = drdy' * B + rxy' * DBdy;  % first order derivatives of shape function to y
end

