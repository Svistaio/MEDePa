function [stimaL,stimaG] = Calcola_Estimatore(ordP,f,phi,dphi,Dphi,dt,nt,istante,datiE)

    global geom

    disp("Calcolo degli estimatori.")

        % Estimatore locale (vettore)
    stimaL = zeros(geom.Nobj.N_ele,1);

        % Coefficienti dell'equazione
    ni    = datiE.ni;
    beta  = datiE.beta;
    sigma = datiE.sigma;


        % Pesi e coordinate dei punti di quadratura
    [xhat, yhat, omega] = Pesi_Nodi_Quad();


    % Calcolo della norma delle funzioni (norm01) e dei gradienti (norm02)
    for e=1:geom.Nobj.N_ele

        [dx,dy] = calcoloDxy(geom.obj.T(e,:));

        areaT = geom.support.TInfo(e).Area;
        
            % Calcolo del lato massimo e delle coordinate dei vertici
        [hE,a1,a2,a3] = LunghezzaLatoMassimo(e);

            % Calcolo della matrice B inversa e inversa trasporta
        Binv  = (1/(2*areaT))*...
                [dy(1) dx(1)
                 dy(2) dx(2)];
    
        BinvT = (1/(2*areaT))*...
                [dy(1) dy(2)
                 dx(1) dx(2)];

            % Funzione FB
        FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

        Re = Calcolo_ResiduoE(omega,xhat,yhat,dt,nt,istante, ...
                              f,FB,phi,dphi,Dphi,Binv,BinvT,...
                              ni,beta,sigma,areaT,ordP,hE,e);

        Se = Calcolo_SaltoE(nt,dphi,BinvT,ni,ordP,e);

            % Stima locale dell'errore sul triangolo e
        stimaL(e) = Re + Se;

    end

        % Estimatore globale (scalare)
    stimaG = sqrt(sum(stimaL));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funzioni ausiliarie %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dx, dy] = calcoloDxy(matrVertT)

    global geom

    % Elementi per la costruzione della matrice B^-T
    dx = zeros(3,1);
    dy = zeros(3,1);

    for i = 1:3
        % Punti i cui nomi rispettano l'ordinamento antiorario
        iS = matrVertT((mod(i+2,3)==0)*3+(mod(i+2,3)~=0)*mod(i+2,3)); % Indice del punto di sinistra
        iD = matrVertT((mod(i+1,3)==0)*3+(mod(i+1,3)~=0)*mod(i+1,3)); % Indice del punto di destra
        pS = geom.obj.P(iS,:); % Punto di sinistra
        pD = geom.obj.P(iD,:); % Punto di destra

        dx(i) = pS(1)-pD(1);
        dy(i) = pD(2)-pS(2);
    end

end

function [hE,a1,a2,a3] = LunghezzaLatoMassimo(e)

    global geom

        % Estrazione coordinate dei tre punti dell'e-esimo triangolo
    a1 = geom.obj.P(geom.obj.T(e,1),:);
    a2 = geom.obj.P(geom.obj.T(e,2),:);
    a3 = geom.obj.P(geom.obj.T(e,3),:);

        % Calcolo della lunghezza massima tra tutti i lati dell'e-esimo triangolo
    hE = max([norm(a1-a2),norm(a2-a3),norm(a3-a1)]);

end

function [Re] = Calcolo_ResiduoE(omega,xhat,yhat,dt,nt,istante, ...
                                 f,FB,phi,dphi,Dphi,Binv,BinvT,...
                                 ni,beta,sigma,areaT,ordP,hE,e)

    global geom u_est

        % Reinizilizzazione della quadratura per il prossimo triangolo
    quadraturaR = 0; % Quadratura del residuo

    for q=1:length(omega)
        Pq = FB(xhat(q),yhat(q)); % Punto di quadratura

        addendoQR = f(Pq(1),Pq(2),nt*dt+istante);
        Phi = phi(xhat(q),yhat(q));
        dPhi = dphi(xhat(q),yhat(q));
        DPhi = reshape(Dphi(xhat(q),yhat(q)),2,2,3*ordP);

        for k=1:3*ordP
            addendoQR =   addendoQR ...
                        + u_est(geom.obj.T(e,k),nt+1)...
                        *(ni*sum((Binv*BinvT).*DPhi(:,:,k),'all') ...
                          -beta'*BinvT*dPhi(:,k) ...
                          -sigma*Phi(k));
        end
        
        quadraturaR = quadraturaR + omega(q)*addendoQR^2;
    end

    Re = (hE^2/ni)*(2*areaT*quadraturaR);

end

function [Se] = Calcolo_SaltoE(nt,dphi,BinvT,ni,ordP,e)

    global geom u_est

    Se = 0;

    RotAnt = [0 -1;1 0]; % Rotazione antioraria = [cos(90°) -sin(90°);sin(90°) cos(90°)]
    RotOra = [0 1;-1 0]; % Rotazione oraria = [cos(90°) sin(90°);-sin(90°) cos(90°)]

        % Ciclo for su ogni lato
    for l = 1:3
        if(geom.obj.Neigh(e,l)>0) % Se il lato non è di bordo
        
                % Reinizilizzazione della quadratura per il prossimo lato
            quadraturaS = 0; % Quadratura del salto

                % Indici del triangolo interno e adiacente
            es = e;
            BinvTs = BinvT;

            ed = geom.obj.Neigh(e,l);
            [dx,dy] = calcoloDxy(geom.obj.T(ed,:));
            areaTd = geom.support.TInfo(ed).Area;
            BinvTd = (1/(2*areaTd))*...
                     [dy(1) dy(2)
                      dx(1) dx(2)];

                % Indice del lato intermedio i due triangoli
            nl = geom.obj.Neigh(e,l+3);
    
            ini = geom.obj.E(nl,1);  % Indice del nodo d'inizio del lato
            inf = geom.obj.E(nl,2);  % Indice del nodo di fine del lato
            
            pi = geom.obj.P(ini,:);  % Coordinate del punto iniziale
            pf = geom.obj.P(inf,:);  % Coordinate del punto finale
            ll = norm(pf-pi,2);      % Lunghezza del lato

            normale = VersoreNorm(ordP,es,ini,inf,pi,pf,RotAnt,RotOra);

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
    
            nodiQ = [0 0.5 1]';
            pesiQ = [1/6 2/3 1/6]';

            dphiPs = ProiettaFunzioneLato(ordP,dphi,es,ini,inf);
            dphiPd = ProiettaFunzioneLato(ordP,dphi,ed,ini,inf);

            qldPhiPs = zeros(6,ordP*3);
            qldPhiPd = zeros(6,ordP*3);
    
            for q=1:3
                qldPhiPs([2*q-1 2*q],:) = dphiPs(nodiQ(q));
                qldPhiPd([2*q-1 2*q],:) = dphiPd(nodiQ(q));
            end
                % S'inverte l'indice associato a «ordP*3» (il secondo) con quello
                % associato a «3» (il terzo): cosí facendo prendendo una faccia
                % del tensore vi sono i gradienti concatenati della medesima funzione
                % di base valutati nei tre punti di quadratura
            qldPhiPs = reshape(qldPhiPs,2,3,ordP*3);
            qldPhiPd = reshape(qldPhiPd,2,3,ordP*3);
            
                % Quadratura del salto
            for q=1:3
                addendoQS = zeros(2,1);
                for k=1:3*ordP
                    addendoQS =   addendoQS ...
                                + ni*u_est(geom.obj.T(es,k),nt+1)*BinvTs*qldPhiPs(:,q,k) ...
                                - ni*u_est(geom.obj.T(ed,k),nt+1)*BinvTd*qldPhiPd(:,q,k);
                end
                quadraturaS = quadraturaS + pesiQ(q)*(dot(addendoQS,normale))^2;        
            end

            Se = Se + (ll/ni)*quadraturaS;
        end
    end

    Se = Se/2;

end

function [n] = VersoreNorm(ordP,nt,ni,nf,pi,pf,RotAnt,RotOra)

    % Funzione per calcolare il versore normale uscente dal triangolo nt;
    % purtroppo è un po' piú complesso di quanto m'aspettavo poiché il
    % triangolatore non ordina in senso antiorario gl'indici globali del
    % lato nella struttura geom.obj.E

    % Dunque bisogna prima accertarsi che la successione dei nodi è
    % antioraria od oraria, ruotando in senso antiorario e orario
    % risepettivamente per poi normalizzare il vettore cosí ottenuto

    global geom

    npT = geom.obj.T(nt,:);      % indici globali del triangolo adiacente
    ili = dot(ni==npT,1:3*ordP); % indice locale iniziale
    ils = mod(ili,3)+1;          % indice locale del nodo successivo [in modo antiorario] a quello iniziale

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

function [fP] = ProiettaFunzioneLato(ordP,f,it,ini,inf)

    % Questa funzione proietta una qualunque funzione f lungo il lato descritto
    % dall'indice iniziale ini e quello finale inf dato l'indice del triangolo it

        [condPerc,condLato] = InfoIntornoLato(ordP,it,ini,inf);

        if(condLato == 1)      % Lato tra il secondo e terzo nodo
            if(condPerc == 2) % Senso di percorrenza orario
                fP = @(t) f(0,t);
            else              % Senso di percorrenza antiorario
                fP = @(t) f(0,1-t);
            end
        elseif(condLato == 2)  % Lato tra il terzo e primo nodo
            if(condPerc == 2) % Senso di percorrenza orario
                fP = @(t) f(1-t,0);
            else              % Senso di percorrenza antiorario
                fP = @(t) f(t,0);
            end
        else                   % Lato tra il primo e secondo nodo
            if(condPerc == 2) % Senso di percorrenza orario
                fP = @(t) f(t,1-t);
            else              % Senso di percorrenza antiorario
                fP = @(t) f(1-t,t);
            end
        end

end

function [condPerc,condLato] = InfoIntornoLato(ordP,nt,ni,nf)

    global geom

    % Questa funzione ricava tre informazioni:
    %   - Un vettore legante indici locali-gloabli
    %   - La condizione legata al senso di percorrenza (1<->antiorario o 2<->orario)
    %   - La condizione legata al lato corrente (1<->x=0/y=t, 2<->y=0/x=t e 3<->y=1-x/x=t)

    legameLG = geom.obj.T(nt,:);      % indici globali del triangolo adiacente
    ili = dot(ni==legameLG,1:3*ordP); % indice locale iniziale
    ils = mod(ili,3)+1;               % indice locale del nodo successivo [in modo antiorario] a quello iniziale
        % Attenzione che «mod(ili,3)» non è «mod(ili,6)» perché si considerano solo i vertici principali

    if(legameLG(ils) == nf)  % Se il nodo successivo è quello finale allora il senso è antiorario
        condPerc = 1;
        condLato = mod(ils,3)+1; % Si prende l'indice del nodo successivo al nodo finale
    else                     % Se il nodo precedente è quello finale allora il senso è orario
        condPerc = 2;
        condLato = ils;      % Si prende l'indice del nodo successivo al primo
    end

end



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