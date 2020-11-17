function [ damages ] = SandyDamagesGivenDeltaGMSL( sandyTable, deltaGMSL )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    tableDeltaGMSL = [0 -0.04 -0.08 -0.10 -0.12 -0.14 -0.16 -0.20 -0.24];
    tableDamageColumns = 8:16;
    
    damages = (interp1(tableDeltaGMSL, sandyTable{:,tableDamageColumns}', deltaGMSL, 'linear', 'extrap'))'; 
 
    
end

