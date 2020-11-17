function [result] = PointsInShp(filename, points, placeName)

    shp = shaperead(filename,'Selector',{@(v1)(strcmp(v1,placeName)),'NAME10'});
    
    
    
    result = inpolygon(points(:,1),points(:,2), shp.X, shp.Y);

end