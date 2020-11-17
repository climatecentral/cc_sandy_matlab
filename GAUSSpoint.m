%goble position to every gauss point of each  integral cell
function  [Xg, Yg,WEI,Ngau]=GAUSSpoint(X1cell,Y1cell,X2cell,Y2cell,Ngauss)
 
%INPUT PARAMETERS
%  Xcell Ycell - coordinations of  rectangle background cell
%   WEI - weight value to each gauss point 
%   Ngau  - the number of gauss point  
% OUTPUT PARAMETERS
%    Xg,Yg - coordinations of Gauss points in each cell
%    Ngauss - exponent number of gauss integration
Xg=[];Yg=[];WEI=[];
if (Ngauss==1)
    WEI=2*2;
    Xg=1/2.*(X1cell+X2cell);
    Yg=1/2.*(Y1cell+Y2cell);
    Ngau=1;
elseif (Ngauss==2)
    a=0.5773502692;
    Xg(1)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(1)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(1)=1;
     Xg(2)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(2)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(2)=1;
     Xg(3)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(3)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(3)=1;
     Xg(4)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(4)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(4)=1;
    Ngau=2*2;
elseif (Ngauss==3) 
    a=0.7745966692;
    W1=0.5555555556;W2=0.8888888889;
     Xg(1)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(1)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(1)=W1*W1;
     Xg(2)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(2)=1/2.*(Y1cell+Y2cell);
    WEI(2)=W1*W2;
     Xg(3)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(3)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(3)=W1*W1;
     Xg(4)=1/2.*(X1cell+X2cell);
    Yg(4)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(4)=W1*W2;
     Xg(5)=1/2.*(X1cell+X2cell);
    Yg(5)=1/2.*(Y1cell+Y2cell);
    WEI(5)=W2*W2;
     Xg(6)=1/2.*(X1cell+X2cell);
    Yg(6)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(6)=W1*W2;
     Xg(7)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(7)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(7)=W1*W1;
     Xg(8)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(8)=1/2.*(Y1cell+Y2cell);
    WEI(8)=W1*W2;
     Xg(9)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(9)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(9)=W1*W1;
     Ngau=3*3;
elseif (Ngauss==4)
    a=0.3399810436;b=0.8611363116;
    W1=0.3478548451;W2=0.6521451549;
    Xg(1)=1/2.*(X1cell+X2cell)-b/2*(X2cell-X1cell);
    Yg(1)=1/2.*(Y1cell+Y2cell)-b/2*(Y2cell-Y1cell);
    WEI(1)=W1*W1;
    Xg(2)=1/2.*(X1cell+X2cell)-b/2*(X2cell-X1cell);
    Yg(2)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(2)=W1*W2;
    Xg(3)=1/2.*(X1cell+X2cell)-b/2*(X2cell-X1cell);
    Yg(3)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(3)=W1*W2;
    Xg(4)=1/2.*(X1cell+X2cell)-b/2*(X2cell-X1cell);
    Yg(4)=1/2.*(Y1cell+Y2cell)+b/2*(Y2cell-Y1cell);
    WEI(4)=W1*W1;
    
     Xg(5)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(5)=1/2.*(Y1cell+Y2cell)-b/2*(Y2cell-Y1cell);
    WEI(5)=W2*W1;
    Xg(6)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(6)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(6)=W2*W2;
    Xg(7)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(7)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(7)=W2*W2;
    Xg(8)=1/2.*(X1cell+X2cell)-a/2*(X2cell-X1cell);
    Yg(8)=1/2.*(Y1cell+Y2cell)+b/2*(Y2cell-Y1cell);
    WEI(8)=W2*W1;
    
    
     Xg(9)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(9)=1/2.*(Y1cell+Y2cell)-b/2*(Y2cell-Y1cell);
    WEI(9)=W1*W1;
    Xg(10)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(10)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(10)=W2*W2;
    Xg(11)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(11)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(11)=W2*W2;
    Xg(12)=1/2.*(X1cell+X2cell)+a/2*(X2cell-X1cell);
    Yg(12)=1/2.*(Y1cell+Y2cell)+b/2*(Y2cell-Y1cell);
    WEI(12)=W2*W1;
    
     Xg(13)=1/2.*(X1cell+X2cell)+b/2*(X2cell-X1cell);
    Yg(13)=1/2.*(Y1cell+Y2cell)-b/2*(Y2cell-Y1cell);
    WEI(13)=W1*W1;
    Xg(14)=1/2.*(X1cell+X2cell)+b/2*(X2cell-X1cell);
    Yg(14)=1/2.*(Y1cell+Y2cell)-a/2*(Y2cell-Y1cell);
    WEI(14)=W1*W2;
    Xg(15)=1/2.*(X1cell+X2cell)+b/2*(X2cell-X1cell);
    Yg(15)=1/2.*(Y1cell+Y2cell)+a/2*(Y2cell-Y1cell);
    WEI(15)=W1*W2;
    Xg(16)=1/2.*(X1cell+X2cell)+b/2*(X2cell-X1cell);
    Yg(16)=1/2.*(Y1cell+Y2cell)+b/2*(Y2cell-Y1cell);
    WEI(16)=W1*W1;
     Ngau=4*4;
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    