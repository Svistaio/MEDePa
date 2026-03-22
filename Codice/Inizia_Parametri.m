function [dominio,sA,sT,sPP] = Inizia_Parametri(datiS,datiD,datiT)

        %%% Inizializzazioni %%%
    areaMax   = datiS.areaMax;
    base      = datiS.base;
    nSim      = datiS.nSim;

    tipoSim   = datiS.tipoSim;
    tipoDomIn = datiD.tipoDom;
    infoDomIn = datiD.infoDom;


        %%% Definizione temporali della simulazione %%%
    sA = []; % Stuttura per la memorizzazione delle aree
    sT = []; % Stuttura per la memorizzazione del tempo

    switch tipoSim
        case 1 % Se si cicla lo spazio
            sA.Max = areaMax*base.^(-(0:nSim-1)');
            sA.Eff = zeros(nSim,1);

            sT.dt = datiT.dt;
            sT.T  = datiT.T;
            sT.nt = floor(sT.T./sT.dt)+(sT.T./sT.dt-floor(sT.T./sT.dt) ~= 0);
            sT.I  = datiT.I;
            sT.E  = datiT.evol;
            
        case 2 % Se si cicla il tempo
            sA.Max = areaMax;
            sA.Eff = 0;

            sT.dt = datiT.dt*base.^(-(0:nSim-1)');
            sT.T  = datiT.T;
            sT.nt = floor(sT.T./sT.dt)+(sT.T./sT.dt-floor(sT.T./sT.dt) ~= 0);
            sT.I  = datiT.I;
            sT.E  = datiT.evol;

    end 
            % Note
        % Il campo «sA.Eff» corrisponde all'area massima effettiva tra tutti gli elementi della triangolazione

        % Il termine «(sT.T/passiT-floor(sT.T/passiT) ~= 0)» aggiunge 1 al numero di passi
        % nel caso in cui si approssima per difetto il con «floor» si aggiunge un passo temporale


        %%% Vettori per il posprocessamento %%%
    sPP = []; % Struttura per la memorizzazione dei dati posprocessorali

    sPP.condA = zeros(nSim,1);
    sPP.dimA  = zeros(nSim,1);
    switch datiS.verifica
        case 1 % Se v'è la verifica si considera l'errore
            sPP.errH1 = zeros(nSim,1);
            sPP.errL2 = zeros(nSim,1);
        case 0 % Se non v'è la verifica si considera l'estimatore
            sPP.stimaG = zeros(nSim,1); % Vettore per le stime globali
            sPP.stimaL = []; % Struttura per le stime locali
    end


        %%% Definizione del dominio %%%
    switch lower(tipoDomIn)
        
        case "punti"
            dominio = infoDomIn;

        case "quadrato"
            dominio = zeros(4,2); % Inizializzazione
    
            % La prima riga è uguale al punto («infoDomIn(2:end)») dato in ingresso
            % e corrisponde all'angolo in basso a sinistra del quadrato
            dominio(1,:) = infoDomIn(2:end);
            dominio(2,:) = infoDomIn(2:end)+[infoDomIn(1) 0]; % Il termine «infoDomIn(1)» è uguale alla lunghezza del quadrato
            dominio(3,:) = infoDomIn(2:end)+[infoDomIn(1) infoDomIn(1)];
            dominio(4,:) = infoDomIn(2:end)+[0 infoDomIn(1)];


        case "rettangolo"
            dominio = zeros(4,2); % Inizializzazione
    
            % La prima riga è uguale al punto («infoDomIn(2:end)») dato in ingresso
            % e corrisponde all'angolo in basso a sinistra del rettangolo
            dominio(1,:) = infoDomIn(3:end);
            dominio(2,:) = infoDomIn(3:end)+[infoDomIn(1) 0]; % Il termine «infoDomIn(1)» è uguale alla lunghezza orizzontale del rettangolo
            dominio(3,:) = infoDomIn(3:end)+[infoDomIn(1) infoDomIn(2)];
            dominio(4,:) = infoDomIn(3:end)+[0 infoDomIn(2)]; % Il termine «infoDomIn(2)» è uguale alla lunghezza verticale del rettangolo


        case "elle"
            dominio = zeros(6,2); % Inizializzazione
    
            % La prima riga è uguale al punto («infoDomIn(2:end)») dato in ingresso
            % e corrisponde all'angolo in basso a sinistra del dominio a L
            dominio(1,:) = infoDomIn(2:end);
            dominio(2,:) = infoDomIn(2:end)+[infoDomIn(1) 0]; % Il termine «infoDomIn(1)» è uguale alla lunghezza del quadrato
            dominio(3,:) = infoDomIn(2:end)+[infoDomIn(1) infoDomIn(1)/2];
            dominio(4,:) = infoDomIn(2:end)+[infoDomIn(1)/2 infoDomIn(1)/2];
            dominio(5,:) = infoDomIn(2:end)+[infoDomIn(1)/2 infoDomIn(1)];
            dominio(6,:) = infoDomIn(2:end)+[0 infoDomIn(1)];
        

        case "cerchio"
            n = infoDomIn(1); % Numero di punti lungo il cerchio
            % Con un po' di trigonometria si ricava che «pi*((n-2)/(2*n))» è l'angolo da sottrarre alla serie «(2*pi/n)*(0:n-1)» per ottenere un poliedro
            % inscritto in un certchio col lato piú basso orizzontale. A parole si procede cosí: si fissa orizzontale l'ultimo lato che descriverà rispetto
            % all'origine un triangolo T con angolo α=2π/n all'origine; l'orizzontalità permette poi d'implicare [con semplice trigonometria] che l'angolo tra il
            % lato piú a destra di T e l'asse orizzontale è proprio π/2-α/2=π/2-π/n=π((n-2)/2n)
            itr = (2*pi/n)*(0:n-1)-pi*((n-2)/(2*n)); % Intervallo di punti
            r = infoDomIn(2); % Raggio del cerchio
            dominio = r*[cos(itr); sin(itr)]'+infoDomIn(3:end); % Punti definitenti il cerchio ove «infoDomIn(3:end)» sono le coordinate del relativo centro


        case "triangolo"
            dominio = zeros(3,2); % Inizializzazione
    
            area   = infoDomIn(1);
            lambda = infoDomIn(2);
    
            % La base e l'altezza del triangolo vengono definiti nel seguente
            % modo: «infoDomIn(1)» corrisponde alla radice dell'area mentre
            % «infoDomIn(2)» al valore regolatore che definisce la base e
            % l'altezza secondo le formule «h=area/lambda» e «base=area*lamda»;
            % cosí procedendo con lambda che tende a zero si ha h divergente
            % e il triangolo diventa una spillo, altrimenti è un triangolo
            % sempre piú schiacciato per valori di lambda crescenti
            
            % L'aspetto chiave è che con tali valori l'area è preservata
            % e uguale a «(h*b)/2=(area/lambda)(area*lamda)/2=(area^2)/2»,
            % ragion per cui l'area si dà in ingresso come un valore
            % sotto radice e potenzialmente moltiplicato per due
    
            b = area*lambda;
            h = area/lambda;
    
            dominio(1,:) = infoDomIn(3:end)+[0 b/2];
            dominio(2,:) = infoDomIn(3:end)-[0 b/2];
            dominio(3,:) = infoDomIn(3:end)+[h 0];


        otherwise
            disp("Forma geometrica non contemplata")

    end

end