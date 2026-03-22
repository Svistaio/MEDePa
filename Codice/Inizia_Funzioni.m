function [sFE,sFV,sFB] = Inizia_Funzioni(datiS,datiE,datiF,datiT)

% (ordP,verifica,evoluzione,ni,beta,sigma,rif,drifx,drify,riff,drift,istante)

        %%% Inizializzazioni %%%
    ni    = datiE.ni;
    beta  = datiE.beta; 
    sigma = datiE.sigma;

        % Forzante e[d eventuale] condizione inziale
    switch datiS.verifica

        case 0 % Se non v'è la verifica si pone uguale alla funzione di riferimento data in ingresso
        
                %%% Soluzione analitica e relativo gradiente %%%
            sFV = []; % Struttura delle funzioni vere
                % La si pone vuota per dare comunque un'uscita alla funzione

                %%% Funzione di Dirichlet e di Neumann %%%
            sFE = []; % Struttura delle funzioni dell'equazione
            sFE.gD = @(x,y,t) datiF.rifu(x,y,t);
            sFE.gN = @(x,y,t) [datiF.drifux(x,y,t)
                               datiF.drifuy(x,y,t)];
    
                %%% Forzante e soluzione iniziale %%%
            sFE.f  = @(x,y,t) datiF.riff(x,y,t);
            sFE.u0 = @(x,y) datiT.drifut(x,y,datiT.I);

        case 1 % Se v'è la verifica la forzante dipende dall'equazione del calore
            
                %%% Soluzione analitica e relativo gradiente %%%
            sFV = []; % Struttura delle funzioni vere
            sFV.uV  = @(x,y,t) datiF.rifu(x,y,t);
            sFV.duV = @(x,y,t) [datiF.drifux(x,y,t)
                                datiF.drifuy(x,y,t)];
        
                %%% Funzione di Dirichlet e di Neumann %%%
            sFE = []; % Struttura delle funzioni dell'equazione
            sFE.gD = @(x,y,t) sFV.uV(x,y,t);
            sFE.gN = @(x,y,t) ni*sFV.duV(x,y,t);
    
            switch datiT.evol

                case 0 % Se non v'è evoluzione, si fissa il tempo e la derivata [∂u/∂t](x,y;istante)=[∂u/∂t](x,y;t̃) compare come correzione della f
                        %%% Forzante e soluzione iniziale %%%
                    sFE.f  = @(x,y,t) -ni*datiF.riff(x,y,t)+beta'*sFV.duV(x,y,t)+sigma*sFV.uV(x,y,t);
                    sFE.u0 = @(x,y) sFV.uV(x,y,datiT.I);

                case 1 % Se v'è evoluzione si considera la derivata ∂u/∂t come contributo della forzante f
                        %%% Forzante e soluzione iniziale %%%
                    sFE.f  = @(x,y,t) datiT.drifut(x,y,t)-ni*datiF.riff(x,y,t)+beta'*sFV.duV(x,y,t)+sigma*sFV.uV(x,y,t);
                    sFE.u0 = @(x,y) sFV.uV(x,y,datiT.I);

            end
    end % Si noti che «riff» nell'ultimo caso è il laplaciano di u

            % Nota
        % Nel controllo dell'evoluzione [nel caso 0] f non dipende da ∂u/∂t poiché si vale la seguente relazione
        % «f-drift(x,y,t̃)  = drift(x,y,t̃)-ni*riff(x,y,t̃)+beta'*du_vera(x,y,t̃)+sigma*u_vera(x,y,t̃)-drift(x,y,t̃)
        %                  = -ni*riff(x,y,t̃)+beta'*du_vera(x,y,t̃)+sigma*u_vera(x,y,t̃)»



        %%% Funzioni di base, relativi gradienti e matrice hessiana %%%
    sFB = []; % Struttura delle funzioni di base

    N1    = @(x,y)  x;
    dN1x  = @(x,y)  1;
    dN1xx = @(x,y)  0;
    dN1xy = @(x,y)  0;
    dN1y  = @(x,y)  0;
    dN1yx = @(x,y)  0;
    dN1yy = @(x,y)  0;

    N2    = @(x,y)  y;
    dN2x  = @(x,y)  0;
    dN2xx = @(x,y)  0;
    dN2xy = @(x,y)  0;
    dN2y  = @(x,y)  1;
    dN2yx = @(x,y)  0;
    dN2yy = @(x,y)  0;

    N3    = @(x,y)  1-x-y;
    dN3x  = @(x,y)  -1;
    dN3xx = @(x,y)  0;
    dN3xy = @(x,y)  0;
    dN3y  = @(x,y)  -1;
    dN3yx = @(x,y)  0;
    dN3yy = @(x,y)  0;

    if(datiS.ordP == 1)
        
        phi1  = @(x,y)  N1(x,y);
        dphi1 = @(x,y) [dN1x(x,y)
                        dN1y(x,y)];
        Dphi1 = @(x,y) [dN1xx(x,y) dN1xy(x,y)
                        dN1yx(x,y) dN1yy(x,y)];

        phi2  = @(x,y)  N2(x,y);
        dphi2 = @(x,y) [dN2x(x,y)
                        dN2y(x,y)];
        Dphi2 = @(x,y) [dN2xx(x,y) dN2xy(x,y)
                        dN2yx(x,y) dN2yy(x,y)];

        phi3  = @(x,y)  N3(x,y);
        dphi3 = @(x,y) [dN3x(x,y)
                        dN3y(x,y)];
        Dphi3 = @(x,y) [dN3xx(x,y) dN3xy(x,y)
                        dN3yx(x,y) dN3yy(x,y)];

        sFB.phi  = @(x,y) [phi1(x,y) phi2(x,y) phi3(x,y)];
        sFB.dphi = @(x,y) [dphi1(x,y) dphi2(x,y) dphi3(x,y)];
        sFB.Dphi = @(x,y) [Dphi1(x,y) Dphi2(x,y) Dphi3(x,y)];

    elseif(datiS.ordP == 2)

        phi1  = @(x,y)  2.*N1(x,y).*(N1(x,y)-1/2);
        dphi1 = @(x,y) [2.*dN1x(x,y).*(N1(x,y)-1/2)+2.*N1(x,y).*dN1x(x,y)
                        2.*dN1y(x,y).*(N1(x,y)-1/2)+2.*N1(x,y).*dN1y(x,y)];
        Dphi1 = @(x,y) [2.*dN1xx(x,y).*(N1(x,y)-1/2)+2.*dN1x(x,y).*dN1x(x,y)+2.*dN1x(x,y).*dN1x(x,y)+2.*N1(x,y).*dN1xx(x,y)...
                        2.*dN1xy(x,y).*(N1(x,y)-1/2)+2.*dN1x(x,y).*dN1y(x,y)+2.*dN1y(x,y).*dN1x(x,y)+2.*N1(x,y).*dN1xy(x,y)
                        2.*dN1yx(x,y).*(N1(x,y)-1/2)+2.*dN1y(x,y).*dN1x(x,y)+2.*dN1x(x,y).*dN1y(x,y)+2.*N1(x,y).*dN1yx(x,y)...
                        2.*dN1yy(x,y).*(N1(x,y)-1/2)+2.*dN1y(x,y).*dN1y(x,y)+2.*dN1y(x,y).*dN1y(x,y)+2.*N1(x,y).*dN1yy(x,y)];

        phi2  = @(x,y)  2.*N2(x,y).*(N2(x,y)-1/2);
        dphi2 = @(x,y) [2.*dN2x(x,y).*(N2(x,y)-1/2)+2.*N2(x,y).*dN2x(x,y)
                        2.*dN2y(x,y).*(N2(x,y)-1/2)+2.*N2(x,y).*dN2y(x,y)];
        Dphi2 = @(x,y) [2.*dN2xx(x,y).*(N2(x,y)-1/2)+2.*dN2x(x,y).*dN2x(x,y)+2.*dN2x(x,y).*dN2x(x,y)+2.*N2(x,y).*dN2xx(x,y)...
                        2.*dN2xy(x,y).*(N2(x,y)-1/2)+2.*dN2x(x,y).*dN2y(x,y)+2.*dN2y(x,y).*dN2x(x,y)+2.*N2(x,y).*dN2xy(x,y)
                        2.*dN2yx(x,y).*(N2(x,y)-1/2)+2.*dN2y(x,y).*dN2x(x,y)+2.*dN2x(x,y).*dN2y(x,y)+2.*N2(x,y).*dN2yx(x,y)...
                        2.*dN2yy(x,y).*(N2(x,y)-1/2)+2.*dN2y(x,y).*dN2y(x,y)+2.*dN2y(x,y).*dN2y(x,y)+2.*N2(x,y).*dN2yy(x,y)];

        phi3  = @(x,y)  2.*N3(x,y).*(N3(x,y)-1/2);
        dphi3 = @(x,y) [2.*dN3x(x,y).*(N3(x,y)-1/2)+2.*N3(x,y).*dN3x(x,y)
                        2.*dN3y(x,y).*(N3(x,y)-1/2)+2.*N3(x,y).*dN3y(x,y)];
        Dphi3 = @(x,y) [2.*dN3xx(x,y).*(N3(x,y)-1/2)+2.*dN3x(x,y).*dN3x(x,y)+2.*dN3x(x,y).*dN3x(x,y)+2.*N3(x,y).*dN3xx(x,y)...
                        2.*dN3xy(x,y).*(N3(x,y)-1/2)+2.*dN3x(x,y).*dN3y(x,y)+2.*dN3y(x,y).*dN3x(x,y)+2.*N3(x,y).*dN3xy(x,y)
                        2.*dN3yx(x,y).*(N3(x,y)-1/2)+2.*dN3y(x,y).*dN3x(x,y)+2.*dN3x(x,y).*dN3y(x,y)+2.*N3(x,y).*dN3yx(x,y)...
                        2.*dN3yy(x,y).*(N3(x,y)-1/2)+2.*dN3y(x,y).*dN3y(x,y)+2.*dN3y(x,y).*dN3y(x,y)+2.*N3(x,y).*dN3yy(x,y)];

        phi4  = @(x,y)  4.*N1(x,y).*N2(x,y);
        dphi4 = @(x,y) [4.*dN1x(x,y).*N2(x,y)+4.*N1(x,y).*dN2x(x,y)
                        4.*dN1y(x,y).*N2(x,y)+4.*N1(x,y).*dN2y(x,y)];
        Dphi4 = @(x,y) [4.*dN1xx(x,y).*N2(x,y)+4.*dN1x(x,y).*dN2x(x,y)+4.*dN1x(x,y).*dN2x(x,y)+4.*N1(x,y).*dN2xx(x,y)...
                        4.*dN1xy(x,y).*N2(x,y)+4.*dN1x(x,y).*dN2y(x,y)+4.*dN1y(x,y).*dN2x(x,y)+4.*N1(x,y).*dN2xy(x,y)
                        4.*dN1yx(x,y).*N2(x,y)+4.*dN1y(x,y).*dN2x(x,y)+4.*dN1x(x,y).*dN2y(x,y)+4.*N1(x,y).*dN2yx(x,y)...
                        4.*dN1yy(x,y).*N2(x,y)+4.*dN1y(x,y).*dN2y(x,y)+4.*dN1y(x,y).*dN2y(x,y)+4.*N1(x,y).*dN2yy(x,y)];

        phi5  = @(x,y)  4.*N3(x,y).*N2(x,y);
        dphi5 = @(x,y) [4.*dN3x(x,y).*N2(x,y)+4.*N3(x,y).*dN2x(x,y)
                        4.*dN3y(x,y).*N2(x,y)+4.*N3(x,y).*dN2y(x,y)];
        Dphi5 = @(x,y) [4.*dN3xx(x,y).*N2(x,y)+4.*dN3x(x,y).*dN2x(x,y)+4.*dN3x(x,y).*dN2x(x,y)+4.*N3(x,y).*dN2xx(x,y)...
                        4.*dN3xy(x,y).*N2(x,y)+4.*dN3x(x,y).*dN2y(x,y)+4.*dN3y(x,y).*dN2x(x,y)+4.*N3(x,y).*dN2xy(x,y)
                        4.*dN3yx(x,y).*N2(x,y)+4.*dN3y(x,y).*dN2x(x,y)+4.*dN3x(x,y).*dN2y(x,y)+4.*N3(x,y).*dN2yx(x,y)...
                        4.*dN3yy(x,y).*N2(x,y)+4.*dN3y(x,y).*dN2y(x,y)+4.*dN3y(x,y).*dN2y(x,y)+4.*N3(x,y).*dN2yy(x,y)];

        phi6  = @(x,y)  4.*N3(x,y).*N1(x,y);
        dphi6 = @(x,y) [4.*dN3x(x,y).*N1(x,y)+4.*N3(x,y).*dN1x(x,y)
                        4.*dN3y(x,y).*N1(x,y)+4.*N3(x,y).*dN1y(x,y)];
        Dphi6 = @(x,y) [4.*dN3xx(x,y).*N1(x,y)+4.*dN3x(x,y).*dN1x(x,y)+4.*dN3x(x,y).*dN1x(x,y)+4.*N3(x,y).*dN1xx(x,y)...
                        4.*dN3xy(x,y).*N1(x,y)+4.*dN3x(x,y).*dN1y(x,y)+4.*dN3y(x,y).*dN1x(x,y)+4.*N3(x,y).*dN1xy(x,y)
                        4.*dN3yx(x,y).*N1(x,y)+4.*dN3y(x,y).*dN1x(x,y)+4.*dN3x(x,y).*dN1y(x,y)+4.*N3(x,y).*dN1yx(x,y)...
                        4.*dN3yy(x,y).*N1(x,y)+4.*dN3y(x,y).*dN1y(x,y)+4.*dN3y(x,y).*dN1y(x,y)+4.*N3(x,y).*dN1yy(x,y)];

        sFB.phi  = @(x,y) [phi1(x,y) phi2(x,y) phi3(x,y) phi4(x,y) phi5(x,y) phi6(x,y)];
        sFB.dphi = @(x,y) [dphi1(x,y) dphi2(x,y) dphi3(x,y) dphi4(x,y) dphi5(x,y) dphi6(x,y)];
        sFB.Dphi = @(x,y) [Dphi1(x,y) Dphi2(x,y) Dphi3(x,y) Dphi4(x,y) Dphi5(x,y) Dphi6(x,y)];

    else
        disp("L'ordine inserito non è contemplato.");
    end

end