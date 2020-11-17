function [w, dwdx, dwdy] = Weight2D(type, para, x,y,xI,yI,dmI)
% EVALUATE WEIGHT FUNCTION
%
% SYNTAX: [w, dwdr, dwdrr] = GaussWeight(type, para, di, dmi)
%
% INPUT PARAMETERS
%    type - Type of weight function
%    para - Weight function parameter
%    x,y   - gauss point coordinates matrix
%    xI,yI  -  nodal point coordinate
%    dmI - Support size
% OUTPUT PARAMETERS
%    w    - Value of weight function at r
%    dwdx - Value of first order derivative of weight function with respect to x at r
%    dwdy - Value of first order derivative of weight function with respect to y at r

r  = sqrt((x-xI).^2+(y-yI).^2)/ dmI;
nnodes_x=length(x);
nnodes_y=length(y);
w=zeros(nnodes_x,nnodes_y);
dwdx=zeros(nnodes_x,nnodes_y);
dwdy=zeros(nnodes_x,nnodes_y);
drdx=zeros(nnodes_x,nnodes_y);
drdy=zeros(nnodes_x,nnodes_y);
dwdr=zeros(nnodes_x,nnodes_y);
dwdr=zeros(nnodes_x,nnodes_y);

for  j=1:nnodes_y
    for i=1:nnodes_x
        if(r(i,j)==0)
            drdx(i,j)=0;
            drdy(i,j)=0;
        else
           drdx(i,j) = x(i,j)/(dmI.^2*r(i,j));
           drdy(i,j)=y(i,j)/(dmI.^2*r(i,j));
        end
        
        
% EVALUATE WEIGHT FUNCTION AND ITS FIRST AND SECOND ORDER OF DERIVATIVES WITH RESPECT r AT r

if (type == 'GAUSS')
   [w(i,j),dwdr(i,j)] = Gauss(para,r(i,j));
elseif (type == 'CUBIC')
   [w(i,j),dwdr(i,j)] = Cubic(r(i,j));
elseif (type == 'SPLI3')
   [w(i,j),dwdr(i,j)] = Spline3(r(i,j));
elseif (type == 'SPLI5')
   [w(i,j),dwdr(i,j)] = Spline5(r(i,j));
elseif (type == 'SPLIB')
   [w(i,j),dwdr(i,j)] = BSpline(dmI/2,r(i,j)); 
elseif (type == 'power')
   [w(i,j),dwdr(i,j)] = power_function(para,r(i,j));   
elseif (type == 'CRBF1')
   [w(i,j),dwdr(i,j)] = CSRBF1(r(i,j));
elseif (type == 'CRBF2')
   [w(i,j),dwdr(i,j)] = CSRBF2(r(i,j));
elseif (type == 'CRBF3')
   [w(i,j),dwdr(i,j)] = CSRBF3(r(i,j));
 elseif (type == 'CRBF4')
   [w(i,j),dwdr(i,j)] = CSRBF4(r(i,j));
 elseif (type == 'CRBF5')
   [w(i,j),dwdr(i,j)] = CSRBF5(r(i,j));
 elseif (type == 'CRBF6')
   [w(i,j),dwdr(i,j)] = CSRBF6(r(i,j));
else
   error('Invalid type of weight function.');
end
    end
end

dwdx  = dwdr .* drdx;
dwdy= dwdr .* drdy;



% shape function
function [w,dwdr] = Gauss(beta,r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
  else
   b2 = beta*beta;
   r2 = r*r;
   eb2 = exp(-b2);

   w     = (exp(-b2*r2) - eb2) / (1.0 - eb2);
   dwdr  = -2*b2*r*exp(-b2*r2) / (1.0 - eb2);
   
end

function [w,dwdr] = Cubic(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
else
   w     = 1-6*r^2+8*r^3-3*r^4;
   dwdr  = -12*r+24*r^2-12*r^3;
  
end

function [w,dwdr] = Spline3(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
elseif (r<=0.5)
   w     = 2/3 - 4*r^2 + 4*r^3;
   dwdr  = -8*r + 12*r^2;
 
else
   w     = 4/3 - 4*r + 4*r^2 - 4*r^3/3;
   dwdr  = -4 + 8*r -4*r^2;
  
end
function [w,dwdr] = Spline5(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 else
   w     = 1-10*r^3+15*r^4-6*r^5;
   dwdr  = -30*r^2 + 60*r^3-30*r^4;
  
end
function [w,dwdr] = BSpline(h,r)%h is the distance between two nodes.
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
elseif (r<=0.5)
   w     = 1/(pi*h^3)*(1-6*r^2+6*r^3);
   dwdr  = 1/(pi*h^3)*(-12*r+18*r^2);
 
else
   w     = 2/(pi*h^3)*(1-r)^3;
   dwdr  = -6/(pi*h^3)*(1-r)^2;
  
end


function [w,dwdr] = power_function(arfa,r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
else
    a2 = arfa*arfa;
   r2 = r*r;
    w     = exp(-r2/a2);
   dwdr  = (-2*r/a2)*exp(-r2/a2);
  
end

function [w,dwdr] = CSRBF2(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^6*(6+36*r+82*r^2+72*r^3+30*r^4+5*r^5);
	dwdr  = 11*r*(r+2)*(5*r^3+15*r^2+18*r+4)*(r-1)^5;

end
function [w,dwdr] = CSRBF1(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^4*(4+16*r+12*r^2+3*r^3);
	dwdr  = -4*(1-r)^3*(4+16*r+12*r^2+3*r^3)+(1-r)^4*(16+24*r+9*r^2);

end
function [w,dwdr] = CSRBF3(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = 1/3+r^2-4/3*r^3+2*r^2*log(r);
	dwdr  = 4*r-4*r^2+4*r*log(r);

end
function [w,dwdr] = CSRBF4(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = 1/15+19/6*r^2-16/3*r^3+3*r^4-16/15*r^5+1/6*r^6+2*r^2*log(r);
	dwdr  = 25/3*r-16*r^2+12*r^3-16/3*r^4+r^5+4*r*log(r);

end

function [w,dwdr] = CSRBF5(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^6*(35*r^2+18*r+3);
	dwdr  =-6*(1-r)^5*(35*r^2+18*r+3)+(1-r)^6*(70*r+18);
   
end

function [w,dwdr] = CSRBF6(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^8*(32*r^3+25*r^2+8*r+1);
	dwdr  =-8*(1-r)^7*(32*r^3+25*r^2+8*r+1)+(1-r)^8*(96*r^2+50*r+8);
   
end