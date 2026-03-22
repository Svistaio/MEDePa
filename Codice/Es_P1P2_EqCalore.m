% Soluzione del problema dell'equazione del calore mediante MEF
% ove si assume di conoscere «a priori» la soluzione esatta

% Imposta la cartella corrente a quella contenente il codice («CreaGrafici.m») qui mostrato
cd(fileparts(matlab.desktop.editor.getActiveFilename));

% Pulizia generale
clear        % Pulisce le variabili
close all    % Chiude tutte le figure aperte
clc          % Pulisce la finestra dei comandi

% global D C M Dd Cd Md
global geom A Ad b ud u_est
global termina


    % Interfaccia per ricavare dall'utente i dati necessari
[datiS,datiE,datiD,datiF,datiT] = Interfaccia_Grafica();

if(termina == 1) % Struttura condizionale per terminare il codice
    close all
    return;
end


if(datiS.verifica == 1) % Se si conosce la soluzione analitica

        % Inizializzazioni dei parametri necessari
    [dominio,sA,sT,sPP] = Inizia_Parametri(datiS,datiD,datiT);

        % Inizializzazione delle funzioni necessarie
    [sFE,sFV,sFB] = Inizia_Funzioni(datiS,datiE,datiF,datiT);

    for a = 1:length(sA.Max)

            % Definizione della triangolazione
        geom = Triangolatore(sA.Max(a),datiS.angoloMin,datiD.lati,dominio);
            
        for t = 1:length(sT.dt)

                % Numero dell'iterazione nella finestra di comando
            disp(append("Svolgimento della ",string(a),"° iterazione in spazio e della ",string(t),"° in tempo."))


                % Calcolo delle matrici relative alla triangolazione e al problema
            [sA.Eff(a),ngdlMax] = Assembla_Matrici(datiS.angoloMin,datiS.ordP, ...
                                                   datiS.stabilita,datiS.massa, ...
                                                   datiE.ni,datiE.beta,datiE.sigma, ...
                                                   sFE.f,sFE.gN,sFB.phi,sFB.dphi,sFB.Dphi, ...
                                                   t,sT.nt(t),sT.dt(t),sT.I);
    
                % Calcolo condizionamento e dimensione della matrice A
            sPP.condA(a) = condest(A);
            sPP.dimA(a) = size(A,1); % Equivale al numero di gdl


                % Evoluzione della soluzione
            Evolutore(datiT.evol,sFE.u0, ...
                      sT.dt(t),sT.nt(t),sT.I, ...
                      ngdlMax,sFE.gD);
    
                % Calcolo dell'errrore
            [sPP.errH1(a+t-1), ...
             sPP.errL2(a+t-1)] = Calcola_Errore(datiS.ordP,sFV.uV,sFV.duV, ...
                                                sFB.phi,sFB.dphi, ...
                                                sT.dt(t),sT.nt(t),sT.I);
                % Da definizione non si possono avere due cicli in tempo e
                % in spazio: se a cicla da 1 fino a 5, t dev'essere 1 e viceversa
                % se quest'ultimo andasse da 1 a 5. Dunque l'espressione
                % «a+t-1» cicla da 1 fino a 5 in ambo i casi precedenti
        
            clc % Pulizia generale della finestra di comando

        end
    end

        % Figura per vedere l'evoluzione della soluzione nell'ultima [migliore] iterazione
    Visualizza_Triangolazione(datiS.ordP,datiS.verifica,size(u_est,1), ...
                              sFV.uV,sT.nt(t),sT.dt(t),sT.I,sT.E)

        % Visualizzazione degli errori
    Grafici_Errore(datiS.ordP,sPP.condA,sPP.dimA, ...
                   sPP.errH1,sPP.errL2,sA.Eff,sT.dt, ...
                   datiS.tipoSim);


else % Se non si conosce la soluzione analitica

        % Inizializzazioni dei parametri necessari
    [dominio,sA,sT,sPP] = Inizia_Parametri(datiS,datiD,datiT);

        % Inizializzazione delle funzioni necessarie
    [sFE,~,sFB] = Inizia_Funzioni(datiS,datiE,datiF,datiT);


    for a = 1:length(sA.Max)

            % Definizione della triangolazione
        geom = Triangolatore(sA.Max(a),datiS.angoloMin,datiD.lati,dominio);
            
        for t = 1:length(sT.dt)

                % Numero dell'iterazione nella finestra di comando
            disp(append("Svolgimento della ",string(a),"° iterazione in spazio e della ",string(t),"° in tempo."))


                % Calcolo delle matrici relative alla triangolazione e al problema
            [sA.Eff(a),ngdlMax] = Assembla_Matrici(datiS.angoloMin,datiS.ordP, ...
                                                   datiS.stabilita,datiS.massa, ...
                                                   datiE.ni,datiE.beta,datiE.sigma, ...
                                                   sFE.f,sFE.gN,sFB.phi,sFB.dphi,sFB.Dphi, ...
                                                   t,sT.nt(t),sT.dt(t),sT.I);
    
                % Calcolo condizionamento e dimensione della matrice A
            sPP.condA(a) = condest(A);
            sPP.dimA(a) = size(A,1); % Equivale al numero di gdl


                % Evoluzione della soluzione
            Evolutore(datiT.evol,sFE.u0, ...
                      sT.dt(t),sT.nt(t),sT.I, ...
                      ngdlMax,sFE.gD);
    
                % Calcolo dell'estimatore
            [sPP.stimaL.(append("S",string(a+t-1))),...
             sPP.stimaG(a+t-1)] = Calcola_Estimatore(datiS.ordP,sFE.f, ...
                                                     sFB.phi,sFB.dphi,sFB.Dphi, ...
                                                     sT.dt(t),sT.nt(t),sT.I,datiE);
                % Da definizione non si possono avere due cicli in tempo e
                % in spazio: se a cicla da 1 fino a 5, t dev'essere 1 e viceversa
                % se quest'ultimo andasse da 1 a 5. Dunque l'espressione
                % «a+t-1» cicla da 1 fino a 5 in ambo i casi precedenti
        
            clc % Pulizia generale della finestra di comando

        end
    end

        % Figura per vedere l'evoluzione della soluzione nell'ultima [migliore] iterazione
    Visualizza_Triangolazione(datiS.ordP,datiS.verifica,size(u_est,1), ...
                              [],sT.nt(t),sT.dt(t),sT.I,sT.E)

        % Visualizzazione degli errori
    Grafici_Estimatore(datiS.ordP,sPP.condA,sPP.dimA, ...
                       sPP.stimaG,sA.Eff,sT.dt, ...
                       datiS.tipoSim);

end

% fprintf('%.2e\n',sPP.errL2);
% disp(" ")
% fprintf('%.2e\n',sPP.errH1);