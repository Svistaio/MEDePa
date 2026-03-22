function [errH1,errL2] = Calcola_Errore(ordP,u_vera,du_vera,phi,dphi,dt,nt,istante)

    global geom u_est

    disp("Calcolo degli errori.")

    % Inizializzazione
    norm01 = 0; % Norma in L2 delle funzioni
    norm02 = 0; % Norma in L2 dei gradienti

    % Pesi e coordinate dei punti di quadratura
    [xhat, yhat, omega] = Pesi_Nodi_Quad();

    % Calcolo della norma delle funzioni (norm01) e dei gradienti (norm02)
    for e=1:geom.Nobj.N_ele

        % Reinizilizzazione della quadratura per il prossimo triangolo
        quadraturaF = 0;
        quadraturaG = 0;

        matrVertT = geom.obj.T(e,:); % Matrice degl'indici dei vertici dell'e-esimo triangolo

        [dx,dy] = calcoloDxy(matrVertT);

        areaT = geom.support.TInfo(e).Area;

        a1 = geom.obj.P(matrVertT(1),:);
        a2 = geom.obj.P(matrVertT(2),:);
        a3 = geom.obj.P(matrVertT(3),:);

        BT = (1/(2*areaT))*...
             [dy(1) dy(2)
              dx(1) dx(2)];

        % Funzione FB
        FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);

        for q=1:length(omega)

            Pq = FB(xhat(q),yhat(q)); % Punto di quadratura

            addendoQF = u_vera(Pq(1),Pq(2),nt*dt+istante);
            Phi = phi(xhat(q),yhat(q));

            addendoQG = du_vera(Pq(1),Pq(2),nt*dt+istante);
            dPhi = dphi(xhat(q),yhat(q));

            for k=1:3*ordP
                addendoQF = addendoQF - u_est(matrVertT(k),nt+1)*Phi(k);
                addendoQG = addendoQG - u_est(matrVertT(k),nt+1)*BT*dPhi(:,k);
            end

            quadraturaF = quadraturaF + omega(q)*addendoQF^2;
            quadraturaG = quadraturaG + omega(q)*dot(addendoQG,addendoQG);

        end

        % Aggiornamento dei valori di norma
        norm01 = norm01 + 2*areaT*quadraturaF;
        norm02 = norm02 + 2*areaT*quadraturaG;
    end

    % Errore commesso in norma H1
    errL2 = sqrt(norm01);
    errH1 = sqrt(norm01+norm02);
    % errH1 = sqrt(norm02);

end

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

%%%%%%%%%%%%%%%%%%%%%
%%% Abbreviazioni %%%
%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%
%%% Codice scartato %%%
%%%%%%%%%%%%%%%%%%%%%%%
%
% Questo ciclo in essenza calcolava la norma in L^2(0,T;Ω) (v. p. 117
% Dispense) per valutare correttamente l'errore non solo in spazio ma anche
% in tempo nel caso in cui ci sia evoluzione; tuttavia applicandolo mi sono
% reso conto che è lievemente maggiore rispetto che al calcolo sopra ove si
% valuta l'errore all'ultimo istante temporale; e siccome il calcolo di
% questa norma è decisamente piú oneroso, per via del ciclo temporale,
% conviene calcolare l'errore come già fatto con minime differenze.
%
% for t = 0:nt % Per ogni istante
%     norm01t = 0;
%     norm02t = 0;
% 
%     % Calcolo della norma delle funzioni (norm01) e dei gradienti (norm02)
%     for e=1:geom.Nobj.N_ele
% 
%         % Reinizilizzazione della quadratura per il prossimo triangolo
%         quadraturaF = 0;
%         quadraturaG = 0;
% 
%         matrVertT = geom.obj.T(e,:); % Matrice degl'indici dei vertici dell'e-esimo triangolo
% 
%         [dx,dy] = calcoloDxy(matrVertT);
% 
%         areaT = geom.support.TInfo(e).Area;
% 
%         a1 = geom.obj.P(matrVertT(1),:);
%         a2 = geom.obj.P(matrVertT(2),:);
%         a3 = geom.obj.P(matrVertT(3),:);
% 
%         BT = (1/(2*areaT))*...
%              [dy(1) dy(2)
%               dx(1) dx(2)];
% 
%         % Funzione FB
%         FB = @(x,y) a1.*x + a2.*y +a3.*(1-x-y);
% 
%         for q=1:length(omega)
% 
%             Pq = FB(xhat(q),yhat(q)); % Punto di quadratura
% 
%             addendoQF = u_vera(Pq(1),Pq(2),t*dt+istante);
%             Phi = phi(xhat(q),yhat(q));
% 
%             addendoQG = du_vera(Pq(1),Pq(2),t*dt+istante);
%             dPhi = dphi(xhat(q),yhat(q));
% 
%             for k=1:3*ordP
%                 addendoQF = addendoQF - u_est(matrVertT(k),t+1)*Phi(k);
%                 addendoQG = addendoQG - u_est(matrVertT(k),t+1)*BT*dPhi(:,k);
%             end
% 
%             quadraturaF = quadraturaF + omega(q)*addendoQF^2;
%             quadraturaG = quadraturaG + omega(q)*dot(addendoQG,addendoQG);
% 
%         end
% 
%         % Aggiornamento dei valori di norma
%         norm01t = norm01t + 2*areaT*quadraturaF;
%         norm02t = norm02t + 2*areaT*quadraturaG;
%     end
%     norm01 = norm01 + dt*norm01t;
%     norm02 = norm02 + dt*norm02t;
% end
% 
% % Errore commesso in norma H1
% errL2 = sqrt(norm01);
% errH1 = sqrt(norm01+norm02);

