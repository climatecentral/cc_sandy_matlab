%Run this script first to build the sandyWarp variable, then run
%"SandyErrorCorrection.m" to 

if ~exist('sandyErrorTable', 'var')
    SandyCorrectNaNElev;
end

minLon =  min(modelLon);
maxLon = max(modelLon);
minLat = min(modelLat);
maxLat = max(modelLat);
%Per Degree
pointsPerDegree = 120;
resX = round(pointsPerDegree * (maxLon-minLon));
resY = round(pointsPerDegree * (maxLat-minLat));

width = maxLon - minLon;
height = maxLat - minLat;

[x,y] = meshgrid(linspace(minLon,maxLon,resX),linspace(minLat,maxLat,resY));



nnodes = size(sandyObservationLat,1);
npoints = resX * resY;


weightType = 'SPLI3';
param=0;



distances = 999999 * ones(resY,resX);
 for i=1:resY
     thisLat = minLat + (i)/resY * height;
     for j=1:resX
         thisLon = minLon + (j)/resX * width;
         closestPoint=[0 0];
         for k=1:nnodes
             sqDist = (sandyObservationLat(k) - thisLat) .^ 2 + (sandyObservationLon(k) - thisLon) .^ 2;
             if sqDist < distances(i,j)
                 distances(i,j) = sqDist;
                 closestPoint = [sandyObservationLon(k) sandyObservationLat(k)];
             end
         end

         distances(i,j) = sqrt(distances(i,j));
     end
 end
 
niter = 50;
radii = zeros(niter,1);
residuals = zeros(niter,size(sandyObservationElev,1));
abs_residuals = zeros(niter,size(sandyObservationElev,1));
sandyWarp = zeros(niter,resY,resX);
 
for iter=6
    
    rad = 0.02 + (iter - 1) * 0.02;
    
    radii(iter) = rad;
    dmI = rad * ones(1, nnodes);


    [PHI, DPHIx, DPHIy] = MLS2DShape(1, nnodes, sandyObservationLon,sandyObservationLat, npoints, x,y, dmI, weightType, rad ); 
   
    sandyWarp(iter,:,:) = reshape(PHI * sandyErrorTable.CC_MODEL_06_MINUS_OBS,resY,resX);  % Approximation function
    
    weightMax = Weight2D(weightType,rad,0,0,0,0,3);
    
    for i=1:size(sandyObservationElev)
        pointX = round((resX - 1) * (sandyObservationLon(i)-minLon) / width  + 1);
        pointY = round((resY - 1) * (sandyObservationLat(i)-minLat) / height  + 1);
      
        pointZ = sandyWarp(iter,pointY,pointX);
        
        residual = sandyErrorTable.MODEL_06_MAXELEV_NAVD88(i) - pointZ - sandyObservationElev(i);
        abs_residual = abs(sandyErrorTable.MODEL_06_MAXELEV_NAVD88(i) - pointZ - sandyObservationElev(i));
        residuals(iter,i) = residual;
        abs_residuals(iter,i) = abs_residual;
    end

    disp([iter,rad,round(mean(residuals(iter,:)),5), round(sqrt(mean(residuals(iter,:).^2)),5)]);
    %mean(abs_residuals(iter,:))
end
            
    

