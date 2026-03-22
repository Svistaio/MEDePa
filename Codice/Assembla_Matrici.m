function [areaMaxEff,ngdlMax] = Assembla_Matrici(angoloMin,ordP,stabilita,massa, ...
                                                 ni,beta,sigma, ...
                                                 f,gN,phi,dphi,DPhi, ...
                                                 t,nt,dt,istante)

    global geom

    % geom.obj.T(e,:) = indici dei nodi dell'e-esimo triangolo
    % geom.obj.P = coordinate di tutt'i nodi della triangolazione

    if(ordP == 1) % Se gli elementi sono P1

        % Inizializzazione dei vettori per la matrice sparsa A
        [iRigA, iColA, valA, contA, ...
         iRigAd,iColAd,valAd,contAd, ...
         ngdlMax,nD] = Inizia_Matrici(angoloMin,stabilita,nt,t);

            % Comunicazione utente
        disp("Calcolo delle matrici e dei vettori.")

        if(stabilita == 0) % Se non s'impone stabilità

                % Costruzione delle matrici A, Ad e del vettore b per i P1 senza stabilità
            AssemblaAP1(ordP,ni,beta,sigma,f,phi,...
                        contA,iRigA,iColA,valA, ...
                        contAd,iRigAd,iColAd,valAd,...
                        ngdlMax,nD,massa, ...
                        t,nt,dt,istante)

        else               % Se s'impone stabilità

                % Costruzione delle matrici A, Ad e del vettore b per i P1 con stabilità
            AssemblaAP1Stab(ordP,ni,beta,sigma,f,phi,dphi,...
                            contA,iRigA,iColA,valA, ...
                            contAd,iRigAd,iColAd,valAd,...
                            ngdlMax,nD, ...
                            t,nt,dt,istante)

        end

            % Comunicazione utente
        disp("Calcolo dei contributi di Neumann per i P1.")

            % Contributo dei lati di Neumann per i P1
        ContributibNP1(ordP,gN,nt,dt,istante)

    elseif(ordP == 2) % Se gli elementi finiti sono P2

            % Si estende la geometria per includere i nodi aggiuntivi
        if(t==1) % Se è la prima iterazione nel tempo
                % Comunicazione utente
            disp("Estensione della geometria per i P2.")
            Estendi_geomP2(); % Si estende la geometria
        end

            % Inizializzazione dei vettori per la matrice sparsa A
        [iRigA, iColA, valA, contA, ...
         iRigAd,iColAd,valAd,contAd, ...
         ngdlMax,nD] = Inizia_Matrici(angoloMin,stabilita,nt,t);

            % Comunicazione utente
        disp("Calcolo delle matrici e dei vettori per i P2.")

        if(stabilita == 0) % Se non s'impone stabilità

                % Costruzione delle matrici A, Ad e del vettore b per i P2 senza stabilità
            AssemblaAP2(ordP,ni,beta,sigma, ...
                        f,phi,dphi,...
                        contA,iRigA,iColA,valA, ...
                        contAd,iRigAd,iColAd,valAd,...
                        ngdlMax,nD,massa, ...
                        t,nt,dt,istante)

        else               % Se s'impone stabilità

                % Costruzione delle matrici A, Ad e del vettore b per i P2 con stabilità
            AssemblaAP2Stab(ordP,ni,beta,sigma, ...
                            f,phi,dphi,DPhi,...
                            contA,iRigA,iColA,valA, ...
                            contAd,iRigAd,iColAd,valAd,...
                            ngdlMax,nD, ...
                            t,nt,dt,istante)
            
        end

            % Comunicazione utente
        disp("Calcolo dei contributi di Neumann.")

            % Contributo dei lati di Neumann per i P2
        ContributibNP2(ordP,gN,nt,dt,istante)

    else
        disp("L'ordine inserito non è contemplato.")
    end

    % Attenzione che il triangolatore fornito dal professore ricorda anche le triangolazioni passate per cui
    % «geom.support.TInfo» contiene anche triangoli di precedenti iterazioni; pertanto per l'area massima
    % bisogna prendere solo le aree del numero di triangoli correnti («geom.Nobj.N_ele»)
    areaMaxEff = max([geom.support.TInfo(1:geom.Nobj.N_ele).Area]); % Funziona anche con «vertcat()»

end


%%% %%% %%% %%% %%% %%% %%% 
%%% Funzioni ausiliari  %%% 
%%% %%% %%% %%% %%% %%% %%% 


    %%% Inizializzazione matrici

function [iRigA,iColA,valA,contA, ...
          iRigAd,iColAd,valAd,contAd, ...
          ngdlMax,nD] = Inizia_Matrici(angoloMin,stabilita,nt,t)

    global geom A D C M S Ad Dd Cd Md Sd b ud

    if(t==1) % Se è la prima iterazione nel tempo
        ngdlSup = (floor(360/angoloMin+1)+1)*geom.Nobj.N_node;     % Limite superiore di gradi di libertà per nodo
        ngdlMax = max(geom.piv.piv); % Numero massimo di gradi di libertà
        b = zeros(ngdlMax,nt+1);
    
        A = spalloc(ngdlMax,ngdlMax,ngdlSup);
        D = spalloc(ngdlMax,ngdlMax,ngdlSup);
        C = spalloc(ngdlMax,ngdlMax,ngdlSup);
        M = spalloc(ngdlMax,ngdlMax,ngdlSup);
    
        nD = length(geom.piv.Di(:,1)); % Numero di nodi di Dirichlet
        Ad = spalloc(ngdlMax,nD,ngdlSup);
        Dd = spalloc(ngdlMax,nD,ngdlSup);
        Cd = spalloc(ngdlMax,nD,ngdlSup);
        Md = spalloc(ngdlMax,nD,ngdlSup);
        ud = zeros(nD,1);
    
        r = 3;
    
        if(stabilita == 1)
            r = r + 1; % Se non v'è stabilità le righe sono tre altrimenti sono quattro
            S  = spalloc(ngdlMax,ngdlMax,ngdlSup);
            Sd = spalloc(ngdlMax,nD,ngdlSup);
        end
    
        iRigA  = zeros(ngdlSup,1); iColA  = zeros(ngdlSup,1); valA  = zeros(ngdlSup,r); contA  = 0;
        iRigAd = zeros(ngdlSup,1); iColAd = zeros(ngdlSup,1); valAd = zeros(ngdlSup,r); contAd = 0;
    else
        ngdlMax = max(geom.piv.piv); % Numero massimo di gradi di libertà
        b = zeros(ngdlMax,nt+1);
        iRigA  = []; iColA  = []; valA  = []; contA  = [];
        iRigAd = []; iColAd = []; valAd = []; contAd = [];
        nD = [];
    end
end


    %%% Funzioni P1

        % Senza stabilità
function [] = AssemblaAP1(ordP,ni,beta,sigma,f,phi, ...
                          contA,iRigA,iColA,valA, ...
                          contAd,iRigAd,iColAd,valAd,...
                          ngdlMax,nD,condMassa, ...
                          t,nt,dt,istante)

    global A D C M Ad Dd Cd Md geom b

        % Pesi e coordinate dei punti di quadratura
    [xhat, yhat, omega] = Pesi_Nodi_Quad();

    if(t==1) % Se è la prima iterazione nel tempo
        for e = 1:geom.Nobj.N_ele % Numero di triangoli
        
                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il calcolo di A, Ad e b."))

                % Calcolo dei (δx1,δx2,δx3), (δy1,δy2,δy3)
            [dx,dy] = CalcoloDxy(e);
        
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo

                % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
                % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

            % Una volta il b coi P1 veniva calcolato approssimando l'integrale
            % considerando f costante sul baricentro del triangolo
            % 
            % Tuttavia ho notato che cosí ragionando s'introduceva un
            % errore costante dell'ordine di 10^-5/10^-6 che è piccolo e di
            % conseguenza rimane invisibile negli ordini di convergenza
            % della maggior parte dei casi di studio, salvo l'ET1, ossia il
            % problema evolutivo in tempo e coi P1 ove non si recuperava
            % l'ordine di convergenza 2 di Crank-Nicolson.
            %
            % Sperimentando ho poi scoperto che anche per le funzioni
            % all'interno dello spazio generato dai P1 compievano un errore
            % superiore all'ordine della precisione di macchina (10^-14-10^-16).
            %
            % La soluzione è semplice: migliorare la quadratura che risolve
            % súbito il problema del piccolo errore commesso; in tal non
            % s'inquinano piú i risultati di tutt'i casi di studio. Questa
            % è la ragione per cui queste due righe di codice sono state eliminate

            % bar = geom.support.TInfo(e).CG;     % Coordinate del baricentro del triangolo e-esimo
            % b(ii,n+1) = b(ii,n+1) + (areaT/3)*f(bar(1),bar(2),n*dt+istante);

            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for j = 1:3*ordP
                        jj = geom.piv.piv(geom.obj.T(e,j));
                        if jj>0 % Se è un grado di libertà
        
                            contA = contA +1;
                            iRigA(contA) = ii;
                            iColA(contA) = jj;
                    
                            valA(contA,1) = (ni/(4*areaT))*(dy(i)*dy(j)+dx(i)*dx(j)); % D
                            valA(contA,2) = (1/6)*(beta(1)*dy(j)+beta(2)*dx(j));      % C
                            valA(contA,3) = (areaT/12)*(1+(i==j));                    % M
    
                        else % Altrimenti è associato a nodo di Dirichlet
        
                            contAd = contAd +1;
                            iRigAd(contAd) = ii;
                            iColAd(contAd) = -jj;
        
                            valAd(contAd,1) = (ni/(4*areaT))*(dy(i)*dy(j)+dx(i)*dx(j)); % Dd
                            valAd(contAd,2) = (1/6)*(beta(1)*dy(j)+beta(2)*dx(j));      % Cd
                            valAd(contAd,3) = (areaT/12)*(1+(i==j));                    % Md
        
                        end % if jj
                    end % for j
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2(i,xhat,yhat,omega,f,areaT,phi,FB,n,dt,istante);
                    end
                end % if ii
            end % for i
        end % for e
    
            % Scartamento degli elementi nulli rimanenti
        interA  = 1:contA;  iRigA  = iRigA(interA);   iColA  = iColA(interA);   valA  = valA(interA,:);
        interAd = 1:contAd; iRigAd = iRigAd(interAd); iColAd = iColAd(interAd); valAd = valAd(interAd,:);
    
            % Definizione delle matrici sparse «D»,«C» e «R»
        D = sparse(iRigA,iColA,valA(:,1),ngdlMax,ngdlMax);
        C = sparse(iRigA,iColA,valA(:,2),ngdlMax,ngdlMax);
        M = sparse(iRigA,iColA,valA(:,3),ngdlMax,ngdlMax);
    
            % Applicazione della stabilità mediante ammassamento
        if(condMassa == 1)
            M = diag(sum(M,2));
        end % Somma delle righe e successiva diagonalizzazione del vettore ottenuto
    
            % Definizione della matrica sparsa «A»
        A  = D + C + sigma*M;
    
            % Definizione delle matrici sparse «Dd»,«Cd» e «Rd»
        Dd = sparse(iRigAd,iColAd,valAd(:,1),ngdlMax,nD);
        Cd = sparse(iRigAd,iColAd,valAd(:,2),ngdlMax,nD);
        Md = sparse(iRigAd,iColAd,valAd(:,3),ngdlMax,nD);
        Ad = Dd + Cd + sigma*Md; % Definizione della matrice sparsa Ad
    else
        for e = 1:geom.Nobj.N_ele % Numero di triangoli

                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il ricalcolo di b."))

            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo

            % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
            % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);
    
            % Una volta il b coi P1 veniva calcolato approssimando l'integrale
            % considerando f costante sul baricentro del triangolo
            % 
            % Tuttavia ho notato che cosí ragionando s'introduceva un
            % errore costante dell'ordine di 10^-5/10^-6 che è piccolo e di
            % conseguenza rimane invisibile negli ordini di convergenza
            % della maggior parte dei casi di studio, salvo l'ET1, ossia il
            % problema evolutivo in tempo e coi P1 ove non si recuperava
            % l'ordine di convergenza 2 di Crank-Nicolson.
            %
            % Sperimentando ho poi scoperto che anche per le funzioni
            % all'interno dello spazio generato dai P1 compievano un errore
            % superiore all'ordine della precisione di macchina (10^-14-10^-16).
            %
            % La soluzione è semplice: migliorare la quadratura che risolve
            % súbito il problema del piccolo errore commesso; in tal non
            % s'inquinano piú i risultati di tutt'i casi di studio. Questa
            % è la ragione per cui queste due righe di codice sono state eliminate

            % bar = geom.support.TInfo(e).CG; % Coordinate del baricentro del triangolo e-esimo
            % b(ii,n+1) = b(ii,n+1) + (areaT/3)*f(bar(1),bar(2),n*dt+istante);
        
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2(i,xhat,yhat,omega,f,areaT,phi,FB,n,dt,istante);
                    end
                end % if ii
            end % for i

        end % for e
    end

end

function [] = ContributibNP1(ordP,gN,nt,dt,istante)

    global geom b

    RotAnt = [0 -1;1 0]; % Rotazione antioraria = [cos(90°) -sin(90°);sin(90°) cos(90°)]
    RotOra = [0 1;-1 0]; % Rotazione oraria = [cos(90°) sin(90°);-sin(90°) cos(90°)]

    for e = 1:length(geom.piv.Ne) % Numero di lati di Neumann

            % Comunicazione utente
        disp(append("Contributo di Neumann del ",string(e),"° lato."))

        nl = geom.piv.Ne(e,1);  % Indice del lato di Neumann considerato
    
        ni = geom.obj.E(nl,1);  % Indice del nodo d'inizio del lato
        nf = geom.obj.E(nl,2);  % Indice del nodo di fine del lato
        
        ji = geom.piv.piv(ni);  % G.d.l. associato al nodo di fine
        jf = geom.piv.piv(nf);  % G.d.l. associato al nodo d'inizio
    
        pi = geom.obj.P(ni,:);       % Coordinate del punto iniziale
        pf = geom.obj.P(nf,:);       % Coordinate del punto finale
    
        normale = VersoreNorm(ordP,nl,ni,nf,pi,pf,RotAnt,RotOra);
    
        l = norm(pf-pi,2);      % Lunghezza del lato
    
        % La formula di quadratura interpolatoria considerata è quella di Simpson
        % (v. p. 280, «Metodi Numerici {27-09-2022}»): [(b-a)/6]*[f(a)+4f((a+b)/2)+f(b)],
        % e come si può vedere essa ragiona in tre nodi di quadratura t1=a, t2=(a+b)/2 e t3=b
        % con pesi w1=(b-a)/6, w2=[2(b-a)]/3 e w3=(b-a)/6; ma in questo contesto a=0 e b=1
        % per cui i nodi diventano t̂1=0, t̂2=1/2 e t̂3=1, invece i pesi sono w1=1/6, w2=2/3 e w3=1/6.
        %
        % Dunque la sommatoria Σ_{q=1}^{Nq}[wq*gN(γe(t̂q))*t̂q] per il nodo d'inizio diventa:
        % Σ_{q=1}^{3}[wq*gN(γe(t̂q))*t̂q]=[(1/6)*gN(pi)*0+(2/3)*gN([pi+pf]/2)*(1/2)+(1/6)*gN(pf)*1]
        %                              =(1/3)*gN([pi+pf]/2)+(1/6)*gN(pf)
        % ove γe(t̂)=pi*(1-t̂)+pf*t̂; d'altro Σ_{q=1}^{Nq}[wq*gN(γe(t̂q))*(1-t̂q)] segue per simmetria:
        % Σ_{q=1}^{3}[wq*gN(γe(t̂q))*(1-t̂q)]=[(1/6)*gN(pi)*1+(2/3)*gN([pi+pf]/2)*(1/2)+(1/6)*gN(pf)*0]
        %                              =(1/6)*gN(pi)+(1/3)*gN([pi+pf]/2)
    
        gammaE = @(t) pi.*(1-t)+pf.*t;
    
        nodiQ = [0 0.5 1]';
        pesiQ = [1/6 2/3 1/6]';
        coordQ = [gammaE(nodiQ(1))
                  gammaE(nodiQ(2))
                  gammaE(nodiQ(3))];
    

        for n = 0:nt

            dNorm = [gN(coordQ(1,1),coordQ(1,2),n*dt+istante)'
                     gN(coordQ(2,1),coordQ(2,2),n*dt+istante)'
                     gN(coordQ(3,1),coordQ(3,2),n*dt+istante)'];
        
            gNQ = dNorm*normale; % Derivata normale imposta perché si conosce «a priori» la soluzione esatta
        
            if(ji>0) % Se il nodo è un g.d.l.
                b(ji,n+1) = b(ji,n+1) + l*sum(pesiQ.*(1-nodiQ).*gNQ); 
            end
            if(jf>0) % Se il nodo è un g.d.l.
                b(jf,n+1) = b(jf,n+1) + l*sum(pesiQ.*nodiQ.*gNQ); 
            end

        end

        % La logica dei contributi di Neumann si basa sul fatto che a ogni
        % g.d.l. viene associata una funzione di base e che nei P1 queste
        % siano funzioni lagrangiane, e quindi di forma nota e semplice.
        %
        % Dalla teoria il contributo del bordo di Neumann si presenta come un
        % integrale di linea da aggiungere al vettore del termine noto b,
        % successivamente scomponibile come somma sui singoli lati; inoltre su
        % questi vivranno piú funzioni di base a seconda dei g.d.l. dei
        % nodi ivi presenti, ma per linearità si possono quindi considerare i
        % contributi delle proriezioni delle singole funzioni di base.
        %
        % Pertanto si prende un lato di Neumann e si legge il nodo iniziale
        % e finale: l'integrale di linea su quel lato segue come coordinate
        % curvilinee quelle che vanno dal nodo iniziale e finale.
        %
        % Ora se sul nodo iniziale v'è un g.d.l. (ji>0) allora si avrà
        % da quadrare una proiezione della relativa funzione di base uguale
        % a 1-t; d'altro canto se v'è sul nodo finale (if>0) si avrà t.

    end

end

        % Con stabilità
function [] = AssemblaAP1Stab(ordP,ni,beta,sigma,f,phi,dphi, ...
                              contA,iRigA,iColA,valA, ...
                              contAd,iRigAd,iColAd,valAd,...
                              ngdlMax,nD, ...
                              t,nt,dt,istante)

    global A D C M S Ad Dd Cd Md Sd geom b

        % Pesi e coordinate dei punti di quadratura
    [xhat, yhat, omega] = Pesi_Nodi_Quad();

    if(t==1) % Se è la prima iterazione nel tempo
        for e = 1:geom.Nobj.N_ele % Numero di triangoli
            
                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il calcolo di A, Ad e b."))

                % Calcolo dei (δx1,δx2,δx3), (δy1,δy2,δy3)
            [dx,dy] = CalcoloDxy(e);
        
                % Calcolo dell'area e delle coordinate baricentriche del triangolo e-esimo
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo

                % Calcolo della matrice B inversa e inversa trasporta
            BinvT = (1/(2*areaT))*...
                    [dy(1) dy(2)
                     dx(1) dx(2)];

                % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
                % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

            % Una volta il b coi P1 veniva calcolato approssimando l'integrale
            % considerando f costante sul baricentro del triangolo
            % 
            % Tuttavia ho notato che cosí ragionando s'introduceva un
            % errore costante dell'ordine di 10^-5/10^-6 che è piccolo e di
            % conseguenza rimane invisibile negli ordini di convergenza
            % della maggior parte dei casi di studio, salvo l'ET1, ossia il
            % problema evolutivo in tempo e coi P1 ove non si recuperava
            % l'ordine di convergenza 2 di Crank-Nicolson.
            %
            % Sperimentando ho poi scoperto che anche per le funzioni
            % all'interno dello spazio generato dai P1 compievano un errore
            % superiore all'ordine della precisione di macchina (10^-14-10^-16).
            %
            % La soluzione è semplice: migliorare la quadratura che risolve
            % súbito il problema del piccolo errore commesso; in tal non
            % s'inquinano piú i risultati di tutt'i casi di studio. Questa
            % è la ragione per cui queste due righe di codice sono state eliminate

            % bar = geom.support.TInfo(e).CG;     % Coordinate del baricentro del triangolo e-esimo
            % b(ii,n+1) = b(ii,n+1) + (areaT/3)*f(bar(1),bar(2),n*dt+istante) ...
            %             + (tauE/2)*(beta(1)*dy(i)+beta(2)*dx(i))*f(bar(1),bar(2),n*dt+istante);


                % Calcolo della costante mk, del lato massimo e del numero di Péclet
            m1  = 1/3;
            hE = LunghezzaLatoMassimo(e);
            PeT = m1*(norm(beta)*hE)/(2*ni);
    
                % Calcolo della viscosità artificiale per stabilizzare la soluzione
            if(PeT<1) % Diffusione dominante
                tauE = m1*hE^2/(4*ni);
            else      % Convezione dominante
                tauE = hE/(2*norm(beta));
            end
    
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for j = 1:3*ordP
                        jj = geom.piv.piv(geom.obj.T(e,j));
                        if jj>0 % Se è un grado di libertà
        
                            contA = contA +1;
                            iRigA(contA) = ii;
                            iColA(contA) = jj;
                    
                            valA(contA,1) = (ni/(4*areaT))*(dy(i)*dy(j)+dx(i)*dx(j));   % D
                            valA(contA,2) = (1/6)*(beta(1)*dy(j)+beta(2)*dx(j));        % C
                            valA(contA,3) = (areaT/12)*(1+(i==j));                      % R
                            valA(contA,4) = (tauE/(4*areaT))*...
                                            (beta(1)*dy(i)+beta(2)*dx(i))*...
                                            (beta(1)*dy(j)+beta(2)*dx(j));              % S
    
                        else % Altrimenti è associato a nodo di Dirichlet
        
                            contAd = contAd + 1;
                            iRigAd(contAd) = ii;
                            iColAd(contAd) = -jj;
        
                            valAd(contAd,1) = (ni/(4*areaT))*(dy(i)*dy(j)+dx(i)*dx(j)); % Dd
                            valAd(contAd,2) = (1/6)*(beta(1)*dy(j)+beta(2)*dx(j));      % Cd
                            valAd(contAd,3) = (areaT/12)*(1+(i==j));                    % Rd
                            valAd(contAd,4) = (tauE/(4*areaT))*...
                                              (beta(1)*dy(i)+beta(2)*dx(i))*...
                                              (beta(1)*dy(j)+beta(2)*dx(j));            % Sd
                            % ud(-jj) = gD(coord(geom.piv.Di(-jj),1),coord(geom.piv.Di(-jj),2));
        
                        end % if jj
                    end % for j
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2Stab(i,xhat,yhat,omega,beta, ...
                                                                  f,areaT,phi,dphi,FB, ...
                                                                  BinvT,tauE,n,dt,istante);                        
                    end
                end % if ii
            end % for i
        end % for e
    
            % Scartamento degli elementi nulli rimanenti
        interA=1:contA;   iRigA=iRigA(interA);    iColA=iColA(interA);    valA=valA(interA,:);
        interAd=1:contAd; iRigAd=iRigAd(interAd); iColAd=iColAd(interAd); valAd=valAd(interAd,:);
    
            % Definizione delle matrici sparse «D»/«Dd»,«C»/«Cd» e «R»/«Rd»
        D = sparse(iRigA,iColA,valA(:,1),ngdlMax,ngdlMax);
        C = sparse(iRigA,iColA,valA(:,2),ngdlMax,ngdlMax);
        M = sparse(iRigA,iColA,valA(:,3),ngdlMax,ngdlMax);
        S = sparse(iRigA,iColA,valA(:,4),ngdlMax,ngdlMax);
    
        % Per mia scelta [momentaneamente] non vi può essere ammassamento assieme alla stabilità
    
            % Definizione della matrice sparsa «A»
        A  = D + C + sigma*M + S;
    
            % Definizione delle matrici sparse «Dd»,«Cd», «Rd» e «Sd»
        Dd = sparse(iRigAd,iColAd,valAd(:,1),ngdlMax,nD);
        Cd = sparse(iRigAd,iColAd,valAd(:,2),ngdlMax,nD);
        Md = sparse(iRigAd,iColAd,valAd(:,3),ngdlMax,nD);
        Sd = sparse(iRigAd,iColAd,valAd(:,4),ngdlMax,nD);
        Ad = Dd + Cd + sigma*Md + Sd; % Definizione della matrice sparsa Ad
    else
        for e = 1:geom.Nobj.N_ele % Numero di triangoli

                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il ricalcolo di b."))

                % Calcolo dei (δx1,δx2,δx3), (δy1,δy2,δy3)
            [dx,dy] = CalcoloDxy(e);

                % Calcolo dell'area e delle coordinate baricentriche del triangolo e-esimo
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo

                % Calcolo della matrice B inversa e inversa trasporta
            BinvT = (1/(2*areaT))*...
                    [dy(1) dy(2)
                     dx(1) dx(2)];

                % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
                % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

            % Una volta il b coi P1 veniva calcolato approssimando l'integrale
            % considerando f costante sul baricentro del triangolo
            % 
            % Tuttavia ho notato che cosí ragionando s'introduceva un
            % errore costante dell'ordine di 10^-5/10^-6 che è piccolo e di
            % conseguenza rimane invisibile negli ordini di convergenza
            % della maggior parte dei casi di studio, salvo l'ET1, ossia il
            % problema evolutivo in tempo e coi P1 ove non si recuperava
            % l'ordine di convergenza 2 di Crank-Nicolson.
            %
            % Sperimentando ho poi scoperto che anche per le funzioni
            % all'interno dello spazio generato dai P1 compievano un errore
            % superiore all'ordine della precisione di macchina (10^-14-10^-16).
            %
            % La soluzione è semplice: migliorare la quadratura che risolve
            % súbito il problema del piccolo errore commesso; in tal non
            % s'inquinano piú i risultati di tutt'i casi di studio. Questa
            % è la ragione per cui queste due righe di codice sono state eliminate

            % bar = geom.support.TInfo(e).CG;     % Coordinate del baricentro del triangolo e-esimo
            % b(ii,n+1) = b(ii,n+1) + (areaT/3)*f(bar(1),bar(2),n*dt+istante) ...
            %             + (tauE/2)*(beta(1)*dy(i)+beta(2)*dx(i))*f(bar(1),bar(2),n*dt+istante);

                % Calcolo della costante mk, del lato massimo e del numero di Péclet
            m1  = 1/3;
            hE = LunghezzaLatoMassimo(e);
            PeT = m1*(norm(beta)*hE)/(2*ni);
    
                % Calcolo della viscosità artificiale per stabilizzare la soluzione
            if(PeT<1) % Diffusione dominante
                tauE = m1*hE^2/(4*ni);
            else      % Convezione dominante
                tauE = hE/(2*norm(beta));
            end
    
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2Stab(i,xhat,yhat,omega,beta, ...
                                                                  f,areaT,phi,dphi,FB, ...
                                                                  BinvT,tauE,n,dt,istante);                        
                    end
                end % if ii
            end % for i
            
        end % for e
    end
end


    %%% Funzioni P2

        % Senza stabilità
function [] = AssemblaAP2(ordP,ni,beta,sigma, ...
                          f,phi,dphi, ...
                          contA,iRigA,iColA,valA, ...
                          contAd,iRigAd,iColAd,valAd,...
                          ngdlMax,nD,condMassa, ...
                          t,nt,dt,istante)

    global A D C M Ad Dd Cd Md geom b

    % Pesi e coordinate dei punti di quadratura
    [xhat, yhat, omega] = Pesi_Nodi_Quad();

    if(t==1) % Se è la prima iterazione nel tempo
        for e = 1:geom.Nobj.N_ele % Numero di triangoli
            
                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il calcolo di A, Ad e b."))

            % Calcolo dei (δx1,δx2,δx3), (δy1,δy2,δy3)
            [dx,dy] = CalcoloDxy(e);
    
                % Calcolo dell'area del triangolo e-esimo
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo
    
            % Calcolo della matrice B inversa e inversa trasporta
            Binv  = (1/(2*areaT))*...
                    [dy(1) dx(1)
                     dy(2) dx(2)];
        
            BinvT = (1/(2*areaT))*...
                    [dy(1) dy(2)
                     dx(1) dx(2)];
    
            % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
            % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);
    
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for j = 1:3*ordP
                        jj = geom.piv.piv(geom.obj.T(e,j));
                        if jj>0 % Se è un grado di libertà
    
                            contA = contA +1;
                            iRigA(contA) = ii;
                            iColA(contA) = jj;
    
                                % Si calcolano in un vettore riga i tre contributi di D, C ed M
                            valA(contA,:) = QuadraturaAP2(i,j,xhat,yhat,omega,ni,beta,phi,dphi,areaT,Binv,BinvT);
    
                        else % Altrimenti è associato a nodo di Dirichlet
    
                            contAd = contAd +1;
                            iRigAd(contAd) = ii;
                            iColAd(contAd) = -jj;
    
                                % Si calcolano in un vettore riga i tre contributi di Dd, Cd ed Md
                            valAd(contAd,:) = QuadraturaAP2(i,j,xhat,yhat,omega,ni,beta,phi,dphi,areaT,Binv,BinvT);
                            % ud(-jj) = gD(coord(geom.piv.Di(-jj),1),coord(geom.piv.Di(-jj),2));
    
                        end % if jj
                    end % for j
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2(i,xhat,yhat,omega,f,areaT,phi,FB,n,dt,istante);
                    end
                end % if ii
            end % for i
        end % for e
    
            % Scartamento degli elementi nulli rimanenti
        interA=1:contA;   iRigA=iRigA(interA);    iColA=iColA(interA);    valA=valA(interA,:);
        interAd=1:contAd; iRigAd=iRigAd(interAd); iColAd=iColAd(interAd); valAd=valAd(interAd,:);
    
            % Definizione delle matrici sparse «D»/«Dd»,«C»/«Cd» e «R»/«Rd»
        D = sparse(iRigA,iColA,valA(:,1),ngdlMax,ngdlMax);
        C = sparse(iRigA,iColA,valA(:,2),ngdlMax,ngdlMax);
        M = sparse(iRigA,iColA,valA(:,3),ngdlMax,ngdlMax);
    
            % Applicazione della stabilità mediante ammassamento
        if(condMassa == 1)
            M = diag(sum(M,2));
        end % Somma delle righe e successiva diagonalizzazione del vettore ottenuto
    
            % Definizione della matrice sparsa «A»
        A  = D + C + sigma*M;
    
            % Definizione delle matrici sparse «Dd»,«Cd» e «Rd»
        Dd = sparse(iRigAd,iColAd,valAd(:,1),ngdlMax,nD);
        Cd = sparse(iRigAd,iColAd,valAd(:,2),ngdlMax,nD);
        Md = sparse(iRigAd,iColAd,valAd(:,3),ngdlMax,nD);
        Ad = Dd + Cd + sigma*Md; % Definizione della matrice sparsa Ad
    else    
        for e = 1:geom.Nobj.N_ele % Numero di triangoli
            
                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il ricalcolo di b."))

                % Calcolo dell'area del triangolo e-esimo
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo
    
            % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
            % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);
    
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2(i,xhat,yhat,omega,f,areaT,phi,FB,n,dt,istante);
                    end
                end % if ii
            end % for i

        end % for e
    end
        
end

function [valA] = QuadraturaAP2(i,j,xhat,yhat,omega,ni,beta,phi,dphi,areaT,Binv,BinvT)

    valA = [0 0 0];

    % % Calcolo della funzione FB se i parametri ni, beta e sigma fossero variabili
    % a1 = geom.obj.P(geom.obj.T(e,1),:);
    % a2 = geom.obj.P(geom.obj.T(e,2),:);
    % a3 = geom.obj.P(geom.obj.T(e,3),:);
    % 
    % % Funzione FB
    % FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

    for q=1:length(omega)

        % Pq = FB(xhat(q),yhat(q)); % Punto di quadratura

        Phi  = phi(xhat(q),yhat(q));
        dPhi = dphi(xhat(q),yhat(q));

        quadratura(1) = ni*dPhi(:,j)'*Binv*BinvT*dPhi(:,i); % D/Dd
        quadratura(2) = beta'*BinvT*dPhi(:,j)*Phi(i);       % C/Cd
        quadratura(3) = Phi(j)*Phi(i);                      % M/Md

        valA = valA + quadratura*omega(q);

    end

    valA = 2*areaT*valA;

end

function [valb] = QuadraturabP2(i,xhat,yhat,omega,f,areaT,phi,FB,n,dt,istante)

    valb = 0;

    for q=1:length(omega)

        Pq = FB(xhat(q),yhat(q)); % Punto di quadratura
        Phi  = phi(xhat(q),yhat(q));

        quadratura = f(Pq(1),Pq(2),n*dt+istante)*Phi(i);
        valb = valb + quadratura*omega(q);

    end

    valb = 2*areaT*valb;

end

function [] = ContributibNP2(ordP,gN,nt,dt,istante)

    global geom b

    RotAnt = [0 -1;1 0]; % Rotazione antioraria = [cos(90°) -sin(90°);sin(90°) cos(90°)]
    RotOra = [0 1;-1 0]; % Rotazione oraria = [cos(90°) sin(90°);-sin(90°) cos(90°)]
    
    for e = 1:length(geom.piv.Ne) % Numero di lati di Neumann

            % Comunicazione utente
        disp(append("Contributo di Neumann del ",string(e),"° lato."))

        nl = geom.piv.Ne(e,1);  % Indice del lato di Neumann considerato
    
        ni = geom.obj.E(nl,1);  % Indice del nodo d'inizio del lato
        nm = geom.obj.E(nl,5);  % Indice del nodo di mezzo del lato [aggiunto alla fine]
        nf = geom.obj.E(nl,2);  % Indice del nodo di fine del lato
        
        ji = geom.piv.piv(ni);  % G.d.l. associato al nodo di fine
        jm = geom.piv.piv(nm);  % G.d.l. associato al nodo di mezzo
        jf = geom.piv.piv(nf);  % G.d.l. associato al nodo d'inizio
    
        pi = geom.obj.P(ni,:);       % Coordinate del punto iniziale
        pf = geom.obj.P(nf,:);       % Coordinate del punto finale
    
        normale = VersoreNorm(ordP,nl,ni,nf,pi,pf,RotAnt,RotOra);
    
        l = norm(pf-pi,2);      % Lunghezza del lato
    
        % La formula di quadratura interpolatoria considerata è quella di Simpson
        % (v. p. 280, «Metodi Numerici {27-09-2022}»): [(b-a)/6]*[f(a)+4f((a+b)/2)+f(b)],
        % e come si può vedere essa ragiona in tre nodi di quadratura t1=a, t2=(a+b)/2 e t3=b
        % con pesi w1=(b-a)/6, w2=[2(b-a)]/3 e w3=(b-a)/6; ma in questo contesto a=0 e b=1
        % per cui i nodi diventano t̂1=0, t̂2=1/2 e t̂3=1, invece i pesi sono w1=1/6, w2=2/3 e w3=1/6.
        %
        % Dunque la sommatoria Σ_{q=1}^{Nq}[wq*gN(γe(t̂q))*t̂q] per il nodo d'inizio diventa:
        % Σ_{q=1}^{3}[wq*gN(γe(t̂q))*t̂q]=[(1/6)*gN(pi)*0+(2/3)*gN([pi+pf]/2)*(1/2)+(1/6)*gN(pf)*1]
        %                              =(1/3)*gN([pi+pf]/2)+(1/6)*gN(pf)
        % ove γe(t̂)=pi*(1-t̂)+pf*t̂; d'altro Σ_{q=1}^{Nq}[wq*gN(γe(t̂q))*(1-t̂q)] segue per simmetria:
        % Σ_{q=1}^{3}[wq*gN(γe(t̂q))*(1-t̂q)]=[(1/6)*gN(pi)*1+(2/3)*gN([pi+pf]/2)*(1/2)+(1/6)*gN(pf)*0]
        %                              =(1/6)*gN(pi)+(1/3)*gN([pi+pf]/2)
    
        gammaE = @(t) pi.*(1-t)+pf.*t;
    
        nodiQ = [0 0.5 1]';
        pesiQ = [1/6 2/3 1/6]';
        coordQ = [gammaE(nodiQ(1))
                  gammaE(nodiQ(2))
                  gammaE(nodiQ(3))];
    
        for n = 0:nt

            dNorm = [gN(coordQ(1,1),coordQ(1,2),n*dt+istante)'
                     gN(coordQ(2,1),coordQ(2,2),n*dt+istante)'
                     gN(coordQ(3,1),coordQ(3,2),n*dt+istante)'];
        
            gNQ = dNorm*normale; % Derivata normale imposta perché si conosce «a priori» la soluzione esatta
        
            if(ji>0) % Se il nodo è un g.d.l.
                b(ji,n+1) = b(ji,n+1) + l*sum(pesiQ.*(2.*(1/2-nodiQ).*(1-nodiQ)).*gNQ); % Funzione di base lagrangiana del nodo iniziale: [(1/2-t)(1-t)]/[(1/2-0)(1-0)]
            end
            if(jm>0) % Se il nodo è un g.d.l.
                b(jm,n+1) = b(jm,n+1) + l*sum(pesiQ.*(4.*nodiQ.*(1-nodiQ)).*gNQ);       % Funzione di base lagrangiana del nodo intermedio: [(0-t)(1-t)]/[(0-1/2)(1-1/2)]
            end
            if(jf>0) % Se il nodo è un g.d.l.
                b(jf,n+1) = b(jf,n+1) + l*sum(pesiQ.*(-2.*nodiQ.*(1/2-nodiQ)).*gNQ);    % Funzione di base lagrangiana del nodo intermedio: [(0-t)(1/2-t)]/[(0-1)(1/2-1)]
            end

        end
        % La logica dei contributi di Neumann si basa sul fatto che a ogni
        % g.d.l. viene associata una funzione di base e che nei P2 queste
        % siano funzioni lagrangiane, e quindi di forma nota e semplice.
        %
        % Dalla teoria il contributo del bordo di Neumann si presenta come un
        % integrale di linea da aggiungere al vettore del termine noto b,
        % successivamente scomponibile come somma sui singoli lati; inoltre su
        % questi vivranno piú funzioni di base a seconda dei g.d.l. dei
        % nodi ivi presenti, ma per linearità si possono quindi considerare i
        % contributi delle proriezioni delle singole funzioni di base.
        %
        % Pertanto si prende un lato di Neumann e si legge il nodo iniziale
        % e finale: l'integrale di linea su quel lato segue come coordinate
        % curvilinee quelle che vanno dal nodo iniziale e finale.
        %
        % Ora se sul nodo iniziale v'è un g.d.l. (ji>0) allora si avrà
        % da quadrare una proiezione della relativa funzione di base uguale a
        % [(1/2-t)(1-t)]/[(1/2-0)(1-0)]; se v'è sul nodo intermedio (jm>0)
        % si avrà [(0-t)(1-t)]/[(0-1/2)(1-1/2)]; infine se v'è sul nodo
        % finale (if>0) si avrà [(0-t)(1/2-t)]/[(0-1)(1/2-1)].

    end

end

        % Con stabilità
function [] = AssemblaAP2Stab(ordP,ni,beta,sigma, ...
                              f,phi,dphi,Dphi, ...
                              contA,iRigA,iColA,valA, ...
                              contAd,iRigAd,iColAd,valAd,...
                              ngdlMax,nD, ...
                              t,nt,dt,istante)

    global A D C M S Ad Dd Cd Md Sd geom b

    % Pesi e coordinate dei punti di quadratura
    [xhat, yhat, omega] = Pesi_Nodi_Quad();
        
    if(t==1) % Se è la prima iterazione nel tempo
        for e = 1:geom.Nobj.N_ele % Numero di triangoli
            
                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo."))

                % Calcolo dei (δx1,δx2,δx3), (δy1,δy2,δy3)
            [dx,dy] = CalcoloDxy(e);
    
                % Calcolo dell'area del triangolo e-esimo
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo
    
                % Calcolo della matrice B inversa e inversa trasporta
            Binv  = (1/(2*areaT))*...
                    [dy(1) dx(1)
                     dy(2) dx(2)];
        
            BinvT = (1/(2*areaT))*...
                    [dy(1) dy(2)
                     dx(1) dx(2)];
    
                % Calcolo della costante mk, del lato massimo e del numero di Péclet
            m2  = 1/24;
            hE = LunghezzaLatoMassimo(e);
            PeT = m2*(norm(beta)*hE)/(2*ni);
    
                % Calcolo della viscosità artificiale per stabilizzare la soluzione
            if(PeT<1) % Diffusione dominante
                tauE = m2*hE^2/(4*ni);
            else      % Convezione dominante
                tauE = hE/(2*norm(beta));
            end
            
                % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
                % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);
    
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for j = 1:3*ordP
                        jj = geom.piv.piv(geom.obj.T(e,j));
                        if jj>0 % Se è un grado di libertà
    
                            contA = contA +1;
                            iRigA(contA) = ii;
                            iColA(contA) = jj;
    
                                % Si calcolano in un vettore riga i quattro contributi di D, C, R ed S
                            valA(contA,:) = QuadraturaAP2Stab(i,j,xhat,yhat,omega, ...
                                                              ni,beta, ...
                                                              phi,dphi,Dphi, ...
                                                              areaT,Binv,BinvT,tauE);
    
                        else % Altrimenti è associato a nodo di Dirichlet
    
                            contAd = contAd +1;
                            iRigAd(contAd) = ii;
                            iColAd(contAd) = -jj;
    
                                % Si calcolano in un vettore riga i quattro contributi di Dd, Cd, Rd ed Sd
                            valAd(contAd,:) = QuadraturaAP2Stab(i,j,xhat,yhat,omega, ...
                                                                ni,beta, ...
                                                                phi,dphi,Dphi, ...
                                                                areaT,Binv,BinvT,tauE);
                            % ud(-jj) = gD(coord(geom.piv.Di(-jj),1),coord(geom.piv.Di(-jj),2));
    
                        end % if jj
                    end % for j
                    for n = 0:nt
                        b(ii) = b(ii) + QuadraturabP2Stab(i,xhat,yhat,omega,beta, ...
                                                          f,areaT,phi,dphi,FB, ...
                                                          BinvT,tauE,n,dt,istante);
                    end
                end % if ii
            end % for i
        end % for e
    
            % Scartamento degli elementi nulli rimanenti
        interA=1:contA;   iRigA=iRigA(interA);    iColA=iColA(interA);    valA=valA(interA,:);
        interAd=1:contAd; iRigAd=iRigAd(interAd); iColAd=iColAd(interAd); valAd=valAd(interAd,:);
    
            % Definizione delle matrici sparse «D»/«Dd»,«C»/«Cd» e «R»/«Rd»
        D = sparse(iRigA,iColA,valA(:,1),ngdlMax,ngdlMax);
        C = sparse(iRigA,iColA,valA(:,2),ngdlMax,ngdlMax);
        M = sparse(iRigA,iColA,valA(:,3),ngdlMax,ngdlMax);
        S = sparse(iRigA,iColA,valA(:,4),ngdlMax,ngdlMax);
    
        % Per mia scelta [momentaneamente] non vi può essere ammassamento assieme alla stabilità
    
            % Definizione della matrice sparsa «A»
        A  = D + C + sigma*M + S;
    
            % Definizione delle matrici sparse «Dd»,«Cd», «Rd» e «Sd»
        Dd = sparse(iRigAd,iColAd,valAd(:,1),ngdlMax,nD);
        Cd = sparse(iRigAd,iColAd,valAd(:,2),ngdlMax,nD);
        Md = sparse(iRigAd,iColAd,valAd(:,3),ngdlMax,nD);
        Sd = sparse(iRigAd,iColAd,valAd(:,4),ngdlMax,nD);
        Ad = Dd + Cd + sigma*Md + Sd; % Definizione della matrice sparsa Ad
    else
        for e = 1:geom.Nobj.N_ele % Numero di triangoli
            
                % Comunicazione utente
            disp(append("Contributo del ",string(e),"° triangolo per il ricalcolo di b."))

                % Calcolo dei (δx1,δx2,δx3), (δy1,δy2,δy3)
            [dx,dy] = CalcoloDxy(e);
    
                % Calcolo dell'area del triangolo e-esimo
            areaT = geom.support.TInfo(e).Area; % Area del triangolo k-esimo
    
                % Calcolo della matrice B inversa e inversa trasporta
            BinvT = (1/(2*areaT))*...
                    [dy(1) dy(2)
                     dx(1) dx(2)];
    
                % Calcolo della costante mk, del lato massimo e del numero di Péclet
            m2  = 1/24;
            hE = LunghezzaLatoMassimo(e);
            PeT = m2*(norm(beta)*hE)/(2*ni);
    
                % Calcolo della viscosità artificiale per stabilizzare la soluzione
            if(PeT<1) % Diffusione dominante
                tauE = m2*hE^2/(4*ni);
            else      % Convezione dominante
                tauE = hE/(2*norm(beta));
            end
            
                % Calcolo della funzione FB se i parametri fossero variabili
            a1 = geom.obj.P(geom.obj.T(e,1),:);
            a2 = geom.obj.P(geom.obj.T(e,2),:);
            a3 = geom.obj.P(geom.obj.T(e,3),:);
        
                % Funzione FB
            FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);
    
            for i = 1:3*ordP
                % Indice cumulativo per la condizione di Dirichlet o Neumann relativo all'indice matrVertT(i) del triangolo
                ii = geom.piv.piv(geom.obj.T(e,i));
                if ii>0   
                    for n = 0:nt
                        b(ii,n+1) = b(ii,n+1) + QuadraturabP2Stab(i,xhat,yhat,omega,beta, ...
                                                                  f,areaT,phi,dphi,FB, ...
                                                                  BinvT,tauE,n,dt,istante);
                    end
                end % if ii
            end % for i

        end % for e
    end

end

function [valA] = QuadraturaAP2Stab(i,j,xhat,yhat,omega, ...
                                    ni,beta, ...
                                    phi,dphi,Dphi, ...
                                    areaT,Binv,BinvT,tauE)

    valA = [0 0 0 0];

    % % Calcolo della funzione FB se i parametri ni, beta e sigma fossero variabili
    % a1 = geom.obj.P(geom.obj.T(e,1),:);
    % a2 = geom.obj.P(geom.obj.T(e,2),:);
    % a3 = geom.obj.P(geom.obj.T(e,3),:);
    % 
    % % Funzione FB
    % FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

    for q=1:length(omega)

        % Pq = FB(xhat(q),yhat(q)); % Punto di quadratura

            % Calcolo delle funzioni di base, dei loro gradienti e delle matrici hessiane nei nodi di quadratura
        Phi  = phi(xhat(q),yhat(q));
        dPhi = dphi(xhat(q),yhat(q));
        DPhi = reshape(Dphi(xhat(q),yhat(q)),2,2,6);
            % La funzione «reshape» riformula l'ultima matrice come collezione di matrici hessiane in un tensore del terz'ordine

        quadratura(1) = ni*dPhi(:,j)'*Binv*BinvT*dPhi(:,i);                     % D/Dd
        quadratura(2) = beta'*BinvT*dPhi(:,j)*Phi(i);                           % C/Cd
        quadratura(3) = Phi(j)*Phi(i);                                          % M/Md
        quadratura(4) = tauE*ni*sum((Binv*BinvT).*DPhi(:,:,j),'all')*(beta'*BinvT*dPhi(:,i)) + ...
                        tauE*(beta'*BinvT*dPhi(:,j))*(beta'*BinvT*dPhi(:,i));   % S/Sd

        valA = valA + quadratura*omega(q);

    end

    valA = 2*areaT*valA;

end

function [valb] = QuadraturabP2Stab(i,xhat,yhat,omega,beta, ...
                                    f,areaT,phi,dphi,FB, ...
                                    BinvT,tauE,n,dt,istante)

    valb = 0;

    for q=1:length(omega)

        Pq = FB(xhat(q),yhat(q)); % Punto di quadratura
        Phi  = phi(xhat(q),yhat(q));
        dPhi = dphi(xhat(q),yhat(q));

        quadratura = f(Pq(1),Pq(2),n*dt+istante)*Phi(i) + ...
                     tauE*f(Pq(1),Pq(2),n*dt+istante)*beta'*BinvT*dPhi(:,i);
        valb = valb + quadratura*omega(q) ;

    end

    valb = 2*areaT*valb;

end



    %%% Funzioni ausiliarie

function [dx,dy] = CalcoloDxy(e)

    global geom

        % Vettori dei vertici degl'indici e delle coordinate dei punti del triangolo
    matrVertT = geom.obj.T(e,:);
    coordP    = geom.obj.P;

        % Elementi per la costruzione della matrice B^-T
    dx = zeros(3,1);
    dy = zeros(3,1);

    for i = 1:3
        % Punti i cui nomi rispettano l'ordinamento antiorario
        iS = matrVertT((mod(i+2,3)==0)*3+(mod(i+2,3)~=0)*mod(i+2,3)); % Indice del punto di sinistra
        iD = matrVertT((mod(i+1,3)==0)*3+(mod(i+1,3)~=0)*mod(i+1,3)); % Indice del punto di destra
        pS = coordP(iS,:); % Punto di sinistra
        pD = coordP(iD,:); % Punto di destra

        dx(i) = pS(1)-pD(1);
        dy(i) = pD(2)-pS(2);
    end

end

function [hE] = LunghezzaLatoMassimo(e)

    global geom

        % Estrazione coordinate dei tre punti dell'e-esimo triangolo
    a1 = geom.obj.P(geom.obj.T(e,1),:);
    a2 = geom.obj.P(geom.obj.T(e,2),:);
    a3 = geom.obj.P(geom.obj.T(e,3),:);

        % Calcolo della lunghezza massima tra tutti i lati dell'e-esimo triangolo
    hE = max([norm(a1-a2),norm(a2-a3),norm(a3-a1)]);

end

function [n] = VersoreNorm(ordP,nl,ni,nf,pi,pf,RotAnt,RotOra)

    % Funzione per calcolare il versore normale uscente dal bordo a seconda
    % del lato considerato; purtroppo è un po' piú complesso di quanto
    % m'aspettavo poiché il triangolatore non ordina in senso antiorario
    % gl'indici globali del lato nella struttura geom.obj.E
    % 
    % Dunque bisogna prima accertarsi che la successione dei nodi è
    % antioraria od oraria, ruotando in senso antiorario e orario
    % risepttivamente per poi normalizzare il vettore cosí ottenuto

    global geom

        % Indice triangolo adiacente al lato di bordo
    nt = geom.obj.E(nl,3:4);
    nt = dot(nt~=-1,nt);

    npT = geom.obj.T(nt,:);      % indici globali del triangolo adiacente
    ili = dot(ni==npT,1:3*ordP); % indice locale iniziale
    ils = mod(ili,3)+1;           % indice locale del nodo successivo [in modo antiorario] a quello iniziale

    if(npT(ils) == nf)       % Se il nodo successivo è quello finale allora il senso è antiorario
        n = RotAnt*(pi-pf)'; % Rotazione antioraria del vettore «pi-pf» orientato in senso orario
    else                     % Se il nodo precedente è quello finale allora il senso è orario
        n = RotOra*(pi-pf)'; % Rotazione oraria del vettore «pi-pf» orientato in senso antiorario
    end
        
        % Ciclo per ingrandire sufficientemente n cosicché la sua normalizzazione non generi problemi
    while(floor(sum(abs(n)))==0)
        n = n*10;
    end
    n = n/norm(n,2);

end


    %%% Funzioni scartate
% function [] = ContributibNStab(ordP,gN,phi,dphi,ni,beta) % Facendo i conti è in realtà inutile
% 
%     Funzione scartata perché l'ipotesi che le funzioni teste venissero
%     modificate nella stabilizzazione della convezione-diffusione è
%     sbagliata in principio: aggiungere la viscosità artificiale non
%     modifica le funzioni teste solo «a posteriori» e limitatamente alla
%     costruzione di A, non nell'integrale di linea per Neumann
%
%     In ogni caso questa funzione è la forma piú generale delle funzioni
%     «ContributibNP1» e «ContributibNP2» perché vale per una generica
%     funzione che poi proietta e integra lungo la corretta linea
% 
%     global geom b
% 
%     RotAnt = [0 -1;1 0]; % Rotazione antioraria = [cos(90°) -sin(90°);sin(90°) cos(90°)]
%     RotOra = [0 1;-1 0]; % Rotazione oraria = [cos(90°) sin(90°);-sin(90°) cos(90°)]
% 
%     for e = 1:length(geom.piv.Ne) % Numero di lati di Neumann
% 
%         il = geom.piv.Ne(e,1);  % Indice del lato di Neumann considerato
% 
%         ini = geom.obj.E(il,1);  % Indice del nodo d'inizio del lato
%         inf = geom.obj.E(il,2);  % Indice del nodo di fine del lato
% 
%         pi = geom.obj.P(ini,:);       % Coordinate del punto iniziale
%         pf = geom.obj.P(inf,:);       % Coordinate del punto finale
% 
%         n = VersoreNorm(ordP,il,ini,inf,pi,pf,RotAnt,RotOra);
% 
%             % Indice triangolo adiacente al lato di bordo
%         it = geom.obj.E(il,3:4);
%         it = dot(it~=-1,it);
% 
%         [legameLG,condPerc,condLato] = InfoIntornoLato(ordP,it,ini,inf);
% 
%         l = norm(pf-pi,2);      % Lunghezza del lato
% 
%         % La formula di quadratura interpolatoria considerata è quella di Simpson
%         % (v. p. 280, «Metodi Numerici {27-09-2022}»): [(b-a)/6]*[f(a)+4f((a+b)/2)+f(b)],
%         % e come si può vedere essa ragiona in tre nodi di quadratura t1=a, t2=(a+b)/2 e t3=b
%         % con pesi w1=(b-a)/6, w2=[2(b-a)]/3 e w3=(b-a)/6; ma in questo contesto a=0 e b=1
%         % per cui i nodi diventano t̂1=0, t̂2=1/2 e t̂3=1, invece i pesi sono w1=1/6, w2=2/3 e w3=1/6.
%         %
%         % Dunque la sommatoria Σ_{q=1}^{Nq}[wq*gN(γe(t̂q))*t̂q] per il nodo d'inizio diventa:
%         % Σ_{q=1}^{3}[wq*gN(γe(t̂q))*t̂q]=[(1/6)*gN(pi)*0+(2/3)*gN([pi+pf]/2)*(1/2)+(1/6)*gN(pf)*1]
%         %                              =(1/3)*gN([pi+pf]/2)+(1/6)*gN(pf)
%         % ove γe(t̂)=pi*(1-t̂)+pf*t̂; d'altro Σ_{q=1}^{Nq}[wq*gN(γe(t̂q))*(1-t̂q)] segue per simmetria:
%         % Σ_{q=1}^{3}[wq*gN(γe(t̂q))*(1-t̂q)]=[(1/6)*gN(pi)*1+(2/3)*gN([pi+pf]/2)*(1/2)+(1/6)*gN(pf)*0]
%         %                              =(1/6)*gN(pi)+(1/3)*gN([pi+pf]/2)
% 
%         gammaE = @(t) pi.*(1-t)+pf.*t;
% 
%         nodiQ = [0 0.5 1]';
%         pesiQ = [1/6 2/3 1/6]';
%         coordQ = [gammaE(nodiQ(1))
%                   gammaE(nodiQ(2))
%                   gammaE(nodiQ(3))];
% 
%         dNorm = [gN(coordQ(1,1),coordQ(1,2))'
%                  gN(coordQ(2,1),coordQ(2,2))'
%                  gN(coordQ(3,1),coordQ(3,2))'];
% 
%         gNQ = dNorm*n; % Derivata normale imposta perché si conosce «a priori» la soluzione esatta
% 
%             % Calcolo della costante mk, del lato massimo e del numero di Péclet
%         if(ordP == 1)
%             m  = 1/3;
%         elseif(ordP == 2)
%             m  = 1/24;
%         else
%             disp("Ordine non contemplato")
%         end
%         hE = LunghezzaLatoMassimo(it);
%         PeT = m*(norm(beta)*hE)/(2*ni);
% 
%             % Calcolo della viscosità artificiale per stabilizzare la soluzione
%         if(PeT<1) % Diffusione dominante
%             tauE = m*hE^2/(4*ni);
%         else      % Convezione dominante
%             tauE = hE/(2*norm(beta));
%         end
% 
% 
%         if(condLato == 1)      % Lato tra il secondo e terzo nodo
%             if(condPerc == 2) % Senso di percorrenza orario
%                 gammaNphi  = @(t) phi(0,t);
%                 gammaNdphi = @(t) dphi(0,t);
%             else              % Senso di percorrenza antiorario
%                 gammaNphi  = @(t) phi(0,1-t);
%                 gammaNdphi = @(t) dphi(0,1-t);
%             end
%         elseif(condLato == 2)  % Lato tra il terzo e primo nodo
%             if(condPerc == 2) % Senso di percorrenza orario
%                 gammaNphi  = @(t) phi(1-t,0);
%                 gammaNdphi = @(t) dphi(1-t,0);
%             else              % Senso di percorrenza antiorario
%                 gammaNphi  = @(t) phi(t,0);
%                 gammaNdphi = @(t) dphi(t,0);
%             end
%         else                   % Lato tra il primo e secondo nodo
%             if(condPerc == 2) % Senso di percorrenza orario
%                 gammaNphi  = @(t) phi(t,1-t);
%                 gammaNdphi = @(t) dphi(t,1-t);
%             else              % Senso di percorrenza antiorario
%                 gammaNphi  = @(t) phi(1-t,t);
%                 gammaNdphi = @(t) dphi(1-t,t);
%             end
%         end
% 
%         gammaNPhi = zeros(3,ordP*3);
%         gammaNdPhi = zeros(6,ordP*3);
% 
%         for k=1:3
%             gammaNPhi(k,:)    = gammaNphi(nodiQ(k));
%             gammaNdPhi([2*k-1 2*k],:) = gammaNdphi(nodiQ(k));
%         end
%             % S'inverte l'indice associato a «ordP*3» (il secondo) con quello
%             % associato a «3» (il terzo): cosí facendo prendendo una faccia
%             % del tensore vi sono i gradienti concatenati della medesima funzione
%             % di base valutati nei tre punti di quadratura
%         gammaNdPhi = reshape(gammaNdPhi,2,3,ordP*3);
%         % [gammaNdphi(nodiQ(1)) ; gammaNdphi(nodiQ(2)) ;  gammaNdphi(nodiQ(3))];
% 
%         for k=1:ordP*3
%             kk = geom.piv.piv(legameLG(k));
%             if(kk>0) % Se il nodo è un g.d.l.
%                 fTeste = gammaNPhi(:,k)+(tauE*beta'*gammaNdPhi(:,:,k))';
%                 b(kk) = b(kk) + l*sum(pesiQ.*fTeste.*gNQ);
%             end
%         end
% 
%         % La logica dei contributi di Neumann si basa sul fatto che a ogni
%         % g.d.l. viene associata una funzione di base e che nei P1 queste
%         % siano funzioni lagrangiane, e quindi di forma nota e semplice.
%         %
%         % Dalla teoria il contributo del bordo di Neumann si presenta come un
%         % integrale di linea da aggiungere al vettore del termine noto b,
%         % successivamente scomponibile come somma sui singoli lati; inoltre su
%         % questi vivranno piú funzioni di base a seconda dei g.d.l. dei
%         % nodi ivi presenti, ma per linearità si possono quindi considerare i
%         % contributi delle proriezioni delle singole funzioni di base.
%         %
%         % Pertanto si prende un lato di Neumann e si legge il nodo iniziale
%         % e finale: l'integrale di linea su quel lato segue come coordinate
%         % curvilinee quelle che vanno dal nodo iniziale e finale.
%         %
%         % Ora se sul nodo iniziale v'è un g.d.l. (ji>0) allora si avrà
%         % da quadrare una proiezione della relativa funzione di base uguale
%         % a 1-t; d'altro canto se v'è sul nodo finale (if>0) si avrà t.
% 
%     end
% 
% end

% function [legameLG,condPerc,condLato] = InfoIntornoLato(ordP,nt,ni,nf)
% 
%     global geom
% 
%     % Questa funzione ricava tre informazioni:
%     %   - Un vettore legante indici locali-gloabli
%     %   - La condizione legata al senso di percorrenza (1<->antiorario o 2<->orario)
%     %   - La condizione legata al lato corrente (1<->x=0/y=t, 2<->y=0/x=t e 3<->y=1-x/x=t)
% 
%     legameLG = geom.obj.T(nt,:);      % indici globali del triangolo adiacente
%     ili = dot(ni==legameLG,1:3*ordP); % indice locale iniziale
%     ils = mod(ili,3)+1;               % indice locale del nodo successivo [in modo antiorario] a quello iniziale
%         % Attenzione che «mod(ili,3)» non è «mod(ili,6)» perché si considerano solo i vertici principali
% 
%     if(legameLG(ils) == nf)  % Se il nodo successivo è quello finale allora il senso è antiorario
%         condPerc = 1;
%         condLato = mod(ils,3)+1; % Si prende l'indice del nodo successivo al nodo finale
%     else                     % Se il nodo precedente è quello finale allora il senso è orario
%         condPerc = 2;
%         condLato = ils;          % Si prende l'indice del nodo successivo al primo
%     end
% 
% end


%%% %%% %%% %%% %%% %%%
%%% Abbreviazioni   %%% 
%%% %%% %%% %%% %%% %%%
% %geom.elements.triangles = TV;
% geom.obj.T = TV;
% 
% %geom.elements.coordinates = V;
% geom.obj.P = V;
% 
% %geom.elements.borders = B;
% geom.obj.E = B;
% 
% %geom.elements.neighbourhood = TT;
% geom.obj.Neigh = TT;
% 
% %geom.elements.vertexesneighbourhood = VB;
% geom.obj.VNeigh = VB;
% 
% %geom.nelements.nTriangles = nT;
% geom.Nobj.N_ele = nT;
% 
% %geom.nelements.nBorders = nB;
% geom.Nobj.N_edge = nB;
% 
% %geom.nelements.nVertexes = nV;
% geom.Nobj.N_node = nV;
% 
% %geom.pivot.nodelist = Nodelist;
% geom.piv.nlist = Nodelist;
% 
% %geom.pivot.Di = Di;
% geom.piv.Di = Di;
% 
% %geom.pivot.Ne = Ne;
% geom.piv.Ne = Ne;
