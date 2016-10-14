classdef DP1
    %DP1 constraint using r-p. dp acts between two markers, from two 
    %separate ridgid bodies
    %   Detailed explanation goes here
    
    properties
        Qi; % [3x1] point location of head of vector AiBar
        Pi; % [3x1] point location of tail of vector AiBar
        Qj; % [3x1] point location of head of vector AjBar
        Pj; % [3x1] point location of tail of vector AjBar
        
    end
    
    %calculated values within
    properties(Dependent)
        AiBar; % vector in local reference frame of body i
        AjBar; % vector in local reference frame of body j
    end
    
    methods
        %constructor
        
        
        %getters fro dependent properties
        function AiBar = get.A
    end
    
end

