classdef Fiber < handle
    %FIBER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        %lunghezza d'onda centrale utilizzata (um)
            lambda = 1.55
        %parametri costruttivi
            length  = 100   %lunghezza della fibra (km)
            n1 = 1.459      %indice di rifrazione del nucleo
            nc = 1.455      %indice di rifrazione del mantello 
            a  = 5          %raggio del nucleo (um)
            b  = 62.5       %raggio del mantello (um)
        
       %parametri di guida
            
            beta2 = -20  % parametro GVD   (ps^2/km)
            beta3 = 0    % parametro 3° ordine (ps^3/km)
            gamma = 0.2  % coeff.  di non linearità (W^-1/km)
            alpha  = 0.2 % coeff. di attenuazione (dB/km) 
    end
    
    methods
        function obj = Fiber(l, beta2, beta3, gamma) 
           if nargin > 0
                obj.length = l;
                obj.beta2 = beta2;
                obj.beta3 = beta3;
                obj.gamma = gamma;
           end
        end
        
       %% _____________________________________________________
       %  Metodi che modificano le proprietà della fibra
        
        function setWavelength (fiber,lambda)
             fiber.lambda= lambda;
        end
        
        function setCoreRadius(fiber,a)
            fiber.a = a;
        end
        
        function setCladdingRadius(fiber, b)
            fiber.b = b;
        end
        
        function setCoreIndex(fiber, n1)
            fiber.n1 = n1;
        end 
        
        function setCladdingIndex(fiber, nc)
            fiber.nc = nc;
        end
        
        function setLength(fiber, l)
            fiber.length = l;
        end
        
        function setBeta2(fiber, beta2)
            fiber.beta2= beta2;
        end
        
        function setBeta3(fiber, beta3)
            fiber.beta3 = beta3;
        end
        
        function setAlpha(fiber, alpha)
            fiber.alpha = alpha;
        end
        
        function setGamma(fiber, gamma)
            fiber.gamma = gamma;
        end
        
        %% ___________________________________________________
        % Metodi che restituiscono le proprietà della fibra
        
        function lambda = getWavelength(fiber)
            lambda = fiber.lambda;
        end
        
        function L = getLength(fiber)
            L = fiber.length;
        end
        
        function a = getCoreRadius(fiber)
            a = fiber.a;
        end
        
        function b = getCladdingRadius(fiber)
            b = fiber.b;
        end
        
        function n1 = getCoreIndex(fiber)
            n1 = fiber.n1;
        end
        
        function nc = getCladdingIndex(fiber)
            nc = fiber.nc;
        end
        
        function beta2 = getBeta2(fiber)
            beta2 = fiber.beta2;
        end 
        
        function beta3 = getBeta3(fiber)
            beta3 = fiber.beta3;
        end
        
        function alpha = getAlpha(fiber)
            alpha = fiber.alpha;
        end
        
        function gamma = getGamma(fiber)
            gamma = fiber.gamma;
        end
        
        function k0 = get_k0(fiber)
            
            k0 = 2*pi/getWavelength(fiber);
        end
        
        function V = get_V(fiber)
            %restituisce il "parametro V"
            k0 = get_k0(fiber);
            n1 = getCoreIndex(fiber);
            nc = getCladdingIndex(fiber);
            a = getCoreRadius(fiber);
            V = k0*a*(n1^2-nc^2)^0.5;
        end
        
        
        function Delta = getDelta(fiber)
            %retituisce la differenza d'indice relativa nucleo-mantello
            n1 = getCoreIndex(fiber);
            nc = getCladdingIndex(fiber);
            Delta = (n1-nc)/n1;
        end
        
        function D = get_D(fiber)
            %restituisce il parametro di dispersione
            beta2 = getBeta2(fiber);
            lambda = getWavelength(fiber);
            k = -1.88365;
            D = k/ lambda^2 * beta2;
        end  
        
        function beta2 = convDtoBeta2(fiber,D)
            %calcola e restituisce beta2 a partire dal parametro D
            lambda = getWavelength(fiber);
            k = -1.88365;
            beta2 = D*lambda^2/k;
        end     
    end
    
    
end
