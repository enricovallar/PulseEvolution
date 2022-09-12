classdef Pulse < handle
    %PULSE realizza degli impulsi che si propagano nelle fibre ottiche 
    %   Detailed explanation goes here
    
    properties
        P0 = 0.1;  %picco di potenza iniziale dell'impulso (W)
        T0 = 30;   %semilarghezza temporale dell'impulso a intensità 1/e (ps) 
        Tmax = 32; %numero massimo di tau
        dtau;      %unità temporale
        tau;       %asse dei tempi normalizzata rispetto a T0
        omega;     %asse delle pulsazioni normalizzate rispetto a T0
        C = 0;     %parametro di chirp in frequenza
        m = 1;     %parametro di pendenza per impulsi super gaussiani
        Uz;        %profilo dell'impulso a distanza z
        type = 'g' %tipo di impulso (g: super-gaussiano, s: a secante iperbolica)
        U0;        %profilo iniziale dell'impulso
        U0_spectrum; %potenza spettrale iniziale dell'impulso 
        max_U0_spec; %massimo valore della potenza spettrale iniziale
        Uz_spectrum; %potenza spettrale a distanza z
        z;           %posizione corrente dell'impulso (km)
        nt = 1024; %punti usati nella FFT
    end
    
    properties(Constant)
        
        
    end
    
    methods
        
        %costruttore 
        function p =  Pulse(C,m,P0, T0,Tmax,  type, nt)
            
            if nargin > 0 
                p.Tmax = Tmax;
            end
            
            if nargin > 6
                p.nt = nt;
            end
            
            
            p.z = 0;                           %l'impulso non e'stato propagato
            p.dtau = (2*p.Tmax)/p.nt;          %asse dei tempi normalizzata rispetto a T0;
            p.tau = [-p.nt/2:p.nt/2-1]*p.dtau; %passo temporale
        
            p.omega = fftshift(-p.nt/2:p.nt/2-1)*(pi/p.Tmax); %asse omega (usato nei calcoli)
            
            %_______________________________________________________________
            %costruttore predefinito
            if nargin == 0                     
                p.U0          = p.generateInputPulse;   %genero il profilo iniziale
                p.Uz = p.U0;                            %inizialmente z=0
                
                %genero lo spettro normalizzato
                temp  = fftshift(ifft(p.Uz));
                p.U0_spectrum = abs(temp).^2;
                p.max_U0_spec = max(p.U0_spectrum);
                p.U0_spectrum = p.U0_spectrum/p.max_U0_spec; %normalizzo;                
                p.Uz_spectrum = p.U0_spectrum; 
                
                return
            end
            %___________________________________________________________
            %costruttore con parametri
            p.C = C;
            p.m = m; 
            p.P0 = P0;
            p.T0 = T0; 
            if nargin > 5
                p.type = type;
            else 
                p.type = 'g'; 
            end  
            
           
            
            p.U0          = p.generateInputPulse; %genero il profilo iniziale
            p.Uz = p.U0;                          %inizialmente z=0
            p.U0_spectrum = ones(1,1024);         %necessario per normalizzare
           
            %genero lo spettro normalizzato
            temp  = fftshift(ifft(p.Uz));
            p.U0_spectrum = abs(temp).^2;
            p.max_U0_spec = max(p.U0_spectrum);
            p.U0_spectrum = p.U0_spectrum/p.max_U0_spec; %normalizzo;                
            p.Uz_spectrum = p.U0_spectrum;   
        end
        
        
        %questa funzione genera il profilo iniziale U0 dell'impulso
        function U0 = generateInputPulse(p)
            if strcmp(p.type, 'g')     %Impulso super-gaussiano 
                U0 = exp(-(1+1i*p.C)*0.5*(p.tau).^(2*p.m));
            elseif strcmp(p.type, 's') %Impulso iperbolico secante
                U0 = sech(p.tau).*exp(-(1i*p.C*p.tau.^2)/2);
            end
        end  
        
        %questa funzione calcola la potenza spettrale dell'impulso a
        %distanza z
        function Uz_spectrum = generateSpectrum(p)
            temp  = fftshift(ifft(p.Uz));
            Uz_spectrum = abs(temp).^2;
            Uz_spectrum = Uz_spectrum/p.max_U0_spec; %normalizzo
        end
        
        %questa funzione restituisce la fase dell'impulso
        function phUz = getPhase(p, uw_c)
            %uw_c : "uw_on"/"uw_off"
            phUz = (angle(p.Uz));
            
            if nargin>1
               if strcmp(uw_c, 'uw_on')
                %elimino i salti aggiungendo multipli di +-2pi
                phUz = unwrap(phUz);   
                phUz = phUz - phUz(p.nt/2+1);%shift rispetto alla fase iniziale
               end   
            end            
        end
        
        %% ____________FUNZIONI DI PLOT____________________________________
        
        
        %rappresento l'impulso sull'asse ax1 e la sua trasformata su ax2
        %con colore c
        function plotPulsePower(p,ax1, ax2,c, l)
            freq = fftshift(p.omega)/(2*pi);
           
            %plot pulse intensity
            if nargin < 5
                plot(ax1, p.tau, abs(p.Uz).^2,'Color',c,'LineStyle','-');
            else
                plot(ax1, p.tau, abs(p.Uz).^2,'Color',c,'LineStyle','-','DisplayName',l);
                legend(ax1,'Location', 'best');
            end 
            ax1.XLimMode = 'auto'; ax1.YLimMode = 'auto';
            ax1.XLabel.String='Normalized Time'; ax1.YLabel.String='Norm. Power';
            ax1.Title.String='Pulse';
            
            %plot pulse spectrum
            if nargin < 5     
                plot(ax2,freq, p.Uz_spectrum, 'Color',c,'LineStyle','-');
            else
                plot(ax2,freq, p.Uz_spectrum, 'Color',c,'LineStyle','-','DisplayName',l);
                legend(ax2,'Location', 'best');
            end   
            ax2.XLimMode = 'auto'; ax2.YLimMode = 'auto';
            ax2.XLabel.String='Normalized Frequency'; ax2.YLabel.String='Spectral Power';
            ax2.Title.String = "Spectrum";
        end
        
        
        %plotta modulo e fase dell'impulso 
        function plotPulseModPhase(p,ax1, ax2,c, l, uw_c)
            
            %mod
            modUz = abs(p.Uz);
            
            %phase 
            phUz = p.getPhase(uw_c);
            
           
            %plot pulse module
            if nargin < 5
                plot(ax1, p.tau, modUz,'Color',c,'LineStyle','-');
            else
                plot(ax1, p.tau, modUz,'Color',c,'LineStyle','-','DisplayName',l);
                legend(ax1,'Location', 'best');
            end
            ax1.XLimMode = 'auto'; ax1.YLimMode = 'auto';
            ax1.XLabel.String='Normalized Time'; ax1.YLabel.String='Normalized Module';
            ax1.Title.String='Pulse Module';
            
            %plot pulse phase
            if nargin < 5     
                plot(ax2,p.tau, phUz, 'Color',c,'LineStyle','-');
            else
                plot(ax2,p.tau, phUz, 'Color',c,'LineStyle','-','DisplayName',l);
                legend(ax2,'Location', 'best');
            end   
            ax2.XLimMode = 'auto'; ax2.YLimMode = 'auto';
            ax2.XLabel.String='Normalized Time'; ax2.YLabel.String='Phase';
            ax2.Title.String = "Pulse Phase";
        end
        
        
        %questo metodo plotta la variazione di frequenza istantanea (freq. chirping)
        function plotPulseFreqChirp(p, ax1,c, l)    
            freqChirp = p.getFreqChirp();
            
            %plot pulse intensity
            if nargin < 4
                plot(ax1, p.tau, freqChirp,'Color',c,'LineStyle','-');
            else
                plot(ax1, p.tau, freqChirp,'Color',c,'LineStyle','-','DisplayName',l);
                legend(ax1,'Location', 'best');
            end
            ax1.XLimMode = 'auto'; ax1.YLimMode = 'auto';
            ax1.XLabel.String='Normalized Time'; ax1.YLabel.String='Frequency Chirp';
            ax1.Title.String='Chirp';
            
        end
 
            
    end

    %%
    %______________________________________________________________________
    %___QUESTI METODI CONTROLLANO L'EVOLUZIONE DELL'IMPULSO
    %______________________________________________________________________
    methods
        %Implementazione dell'algoritmo split step fourier method
        %permette il calcolo della dispersione al terzo ordine 
        function [Uz,powUz_spectrum] =  SSFM(p,f,L,c)
            %>> Input:
            %- p : impulso propagato
            %- f : fibra simulata
            %- L : distanza di propagazione (km)
            %- c : comando "TOD_on"/"TOD_off"
            %      per abilitare la simulazione della dispersione 
            %      con effetti del terzo ordine
            %>> Output:
            %- Uz         : Nuova forma dell'impulso dopo la propagazione
            %- powUz_spec : Potenza spettrale dopo la propagazione 
            
            if(L>f.length)
               error('Hai superato la lunghezza della fibra');
            end
            
            
            %definisco il parametro LD
            L_D  = getL_D(p,f);  %L_D  Lunghezza di dispersione 
            L_DT = getL_DT(p,f); %L_D' Lunghezza di dispersione
                                 %     dovuta a beta3 (TOD par.) 
            L_NL = getL_NL(p,f); %L_NL Lunghezza di nonlinearita'
            
            %definisco il parametro N
            N = getN(p,f);    %N^2 = L_D/L_NL 
            NT= getNT(p,f);   %N^2 = L_D'/L_NL da usare se beta2 = 0;
            
            %ottengo l'attenuazione in (Np/km)
            alpha_Np = f.alpha/4.34;
            
            %definisco il numero di passi 
            if L_NL < 1 
                step_num = round(40*L/L_NL);
            else
                step_num = round(40*L);
            end
            
            %limite inferiore
            if step_num < 40
               step_num = 40;
            end
                 
            %limite superiore
            if step_num > 20000
                step_num = 20000;
            end
            
            %calcolo la lunghezza di un passo (km)
            deltaz = L/step_num;
            
            
            %______________________________________________________________
            %___________dispersione________________________________________
           
            %effetto di beta2 (Group Velocity Dispersion)
            GVD = 0.5/L_D*1i*sign(f.beta2)*p.omega.^2*deltaz; 
            %effetto di beta3 (Third Order Dispersione)
            TOD = 1/(6*L_DT)*1i*sign(f.beta3)*p.omega.^3*deltaz;
            %attenuazione
            ATT = alpha_Np/2*deltaz;
            
            if strcmp(c, "TOD_on")
                %abilito gli effetti di TOD
                dispersion = exp(GVD + TOD - ATT);
            else
                %disabilito gli effetti di TOD
                dispersion = exp(GVD - ATT);
            end 
            
            %______________________________________________________________
            %____________nonlinearita'_____________________________________
            if f.beta2==0 && f.beta3 == 0
                %Caso limite, N = inf
                hhz = 1i/L_NL*deltaz;
            elseif f.beta2 == 0  
                %solo TOD : si usano i parametri N' e L_D'
                hhz = 1i*NT^2/L_DT*deltaz;
            else
                %si usano i parametri N e L_D
                hhz = 1i*N^2/L_D*deltaz;
            end
       
            %______________________________________________________________
            %_____________LOOP PRINCIPALE__________________________________
            
            % Schema:  1/2 N -> D -> 1/2 N
            
            temp = p.U0.*exp(-abs(p.U0).^2.*hhz/2); %1/2 passo iniziale
            
            for n=1 : step_num                      %loop
                f_temp = ifft(temp).*dispersion;    %dispersione
                p.Uz = fft(f_temp);
                temp=p.Uz.*exp(abs(p.Uz).^2.*hhz);  %nonlinearita'
            end
            
            %Risultato 
            p.Uz = temp.*exp(abs(p.Uz).^2.*hhz/2);  %1/2 passo finale
            p.Uz = p.Uz/exp(-alpha_Np*L/2);         %Normalizzo considerando l'attenuazione
            p.Uz_spectrum = p.generateSpectrum();   %potenza spettrale 
            p.z = L;                                %distanza raggiunta
            %output
            Uz = p.Uz;                              %ampiezza normalizzata
            powUz_spectrum = p.Uz_spectrum;         %potenza spettrale
        end 
        
        
        
        
        %%
        %questo metodo calcola la lunghezza di dispersione L_D
        function L_D  = getL_D(p,f)
            L_D = p.T0^2/(abs(f.beta2));
        end
        
        %questo metodo calcola la lunghezza di dispersione L_D
        function L_DT  = getL_DT(p,f)
            L_DT = p.T0^3/(abs(f.beta3));
        end
        
        %questo metodo calcola la lunghezza di nonlinearità L_NL
        function L_NL = getL_NL(p,f)
            L_NL = 1/(p.P0*f.gamma);
        end
        
        %questo metodo calcola il parametro N 
        function N = getN(p,f)
            N = sqrt(f.gamma*p.P0*p.T0^2/abs(f.beta2));
        end
        
        %questo metodo calcola il parametro NT
        function NT = getNT(p,f)
            NT = sqrt(f.gamma*p.P0*p.T0^3/abs(f.beta3));
        end
        
        
        %questo metodo calcola l'evoluzione del parametro di chirp
        function freqChirp = getFreqChirp(p)
            phUz = unwrap(p.getPhase());
            freqChirp=zeros(1,p.nt);
            for i = 1 : (length(phUz)-1)
                freqChirp(i) = -(phUz(i+1)-phUz(i))/p.dtau;
            end
            
        end
        
        %%_________________________________________________________________
        %_____________FULL EVOLUTION_______________________________________
        
        %I metodi seguenti sono utili ad analizzare l'evoluzione
        %dell'impulso lungo tutta la lunghezza della fibra
        
        %questa funzione permette restituisce 
        %la potenza dell'impulso p a diverse
        %distanze tra l'inizio e la fine del tratto di fibra f
        %il numero di campioni successivi all'impulso di input è step_num
        function [z,powUz, powUz_spec] = evolve(p,f,c, step_num)
            deltaL  = f.length/(step_num);
            z = [0:deltaL:f.length];
            powUz      = zeros(step_num+1,p.nt);
            powUz_spec = zeros(step_num+1,p.nt);
            powUz(1, :)      = abs(p.U0).^2;
            powUz_spec(1, :) = abs(p.U0_spectrum).^2;
            for  n = 2 : step_num+1
                [Uz(n,:),powUz_spec(n,:)] = SSFM(p,f,z(n),c); 
                powUz(n,:) = abs(Uz(n,:)).^2;
            end
            
        end
        
        %questo metodo calcola il fattore di allargamento temporale
        function [z,bFactor, s_bFactor] = getBFactor(p,f,c, step_num)
            %evolvo l'impulso lungo la fibra
            [z,powUz, powUz_spec] = evolve(p,f,c,step_num);
            
            %______________________________________________________________
            %_________ALLARGAMENTO TEMPORALE_______________________________
            %larghezza iniziale
            RMS_width_0 = sqrt(var(p.tau*p.T0,powUz(1,:)));
            
            %calcolo del fattore di allargamento
            bFactor = zeros(1,step_num+1);
            for  n = 1:step_num+1
                if  any(isnan(powUz(n,:)),'all')
                    bFactor(n) = NaN;
                else
                    RMS_width_z = sqrt(var(p.tau*p.T0,powUz(n,:)));
                    bFactor(n) = RMS_width_z/RMS_width_0;
                end
            end
            
            %______________________________________________________________
            %________ALLARGAMENTO SPETTRALE________________________________
            %larghezza spettrale iniziale
            freq = fftshift(p.omega)/(2*pi); %asse frequenza
            %compensazione errori numerici:
            [~,powUz_spec_a0] = SSFM(p,f,0,c); %propagazione vicina all'origine
            s_RMS_width_0 = sqrt(var(freq,powUz_spec_a0));
            
            %calcolo del fattore di allargamento
            s_bFactor = zeros(1,step_num+1);
            s_bFactor(1) = 1;
            for  n = 2:step_num+1
                if  any(isnan(powUz_spec(n,:)),'all')
                    s_bFactor(n) = NaN;
                else
                    s_RMS_width_z = sqrt(var(freq, powUz_spec(n,:)));
                    s_bFactor(n) = s_RMS_width_z/s_RMS_width_0;
                end

            end 
        end

        
        
    end
    
end

