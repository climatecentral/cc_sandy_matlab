sandyElevationData = csvread('sandy/maxele_m7_Sandy_ClimateCentral.csv',1,0);
disp('csv done reading');
elevDataLon = sandyElevationData(:,2);
elevDataLat = sandyElevationData(:,3);
elevDataDepth = -1 * sandyElevationData(:,4);

sandyElev = zeros(size(sandyZ));

for i=1:size(sandyZ,1)
    disp(i)
    closestDistance = 99999;
    closestElev = -1;
    
    for j=1:size(elevDataDepth,1)
        latDist = sandyLat(i) - elevDataLat(j);
        lonDist = sandyLon(i) - elevDataLon(j);
        
        dist = norm([latDist,lonDist],2);
        
        if dist < closestDistance
            closestDistance = dist;
            closestElev = elevDataDepth(j);
        end
    end
        
    sandyElev(i) = closestElev;
end