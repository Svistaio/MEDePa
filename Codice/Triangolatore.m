function [geom] = Triangolatore(areaMax,angoloMin,lati,dominio)

    global geom

    addpath('../../Esercizi/Triangolatore/bbtr30_short')
    % disp('../../Esercizi/Triangolatore/bbtr30_short added to the path')
    
    %----------------------------------------------------------------------------
    %
    % Triangolazione di un dominio quadrata
    % con condizioni di Dirichlet sul bordo
    %
    %----------------------------------------------------------------------------
    %
    %  Autore: Stefano Berrone
    %  Politecnico di Torino
    %
    %----------------------------------------------------------------------------
    
    %clc
    %clear all
    
    % -------------------------------
    % Inserimento dei vertici
    % -------------------------------
    
    Domain.InputVertex = dominio;
    
    
    % ---------------------------------------------
    % Definizione del dominio a partire dai Vertici
    % ---------------------------------------------
    
    % Dichiaro le variabili per delimitare il dominio
    Domain.Boundary.Values = 1:length(dominio);
    % lato di bordo 1 dal nodo 1 al nodo 2
    % lato di bordo 2 dal nodo 2 al nodo 3
    % lato di bordo 3 dal nodo 3 al nodo 4
    % lato di bordo 4 dal nodo 4 al nodo 1
    
    Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
    Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio
    
    % --------------------------------------------------
    % Definizione delle condizioni al contorno a partire
    % dai Vertici e dai lati di bordo
    % --------------------------------------------------
    
    % valori numerici per le condizioni al contorno
    BC.Values = [0.0 12.0 0.0 14.0 0.0 16.0 0.0 0.0 0.0];
    % Ho notato che questi sono del tutto irrilevanti nella soluzione
    
    % Se i marcatori sono uguali si riferiscono alla medesima condizione altrimenti si distinguono
    % Ad esempio qui ai vertici si associa la stessa condizione di Neumann
    
    % Di solito Dirichlet s'impone ai vertici mentre Neumann ai lati: infatti
    % la prima è forte e vincola la geometria mentre la seconda è piú larga.
    % (Bordi di Dirichlet chiusi mentre quelli di Neumann sono aperti)
    
    % Il marcatore zero dev'essere messo tra due lati associati a Neumann (apparentemente sempre)
    % I punti interni invece presentano sempre marcatore nullo.
    
    % Marcatori delle condizioni al contorno sui bordi del dominio: dispari -> Dirichlet; pari -> Neumann

        % Lati
    BC.Boundary.Values = lati;

        % Vertici
    vertici = zeros(size(dominio,1),1);
    l = length(lati); 
    for i=1:l
        if(mod(lati(i),2) == 1)     % Lato di Dirichlet
            vertici(i) = 1;
            vertici(mod(i,l)+1) = 1;
        else                        % Lato di Neumann
            
            % Si controllano i due punti del lato: se sono già di Dirichlet si lasciano invariati, altrimenti s'impone nullo,
            % ossia non si dà alcuna condizione perché, come afferma perfettamente il manuale BBTR (p. 12):
            % «These markers must be odd numbers, corresponding to Dirichlet boundary conditions, or zeros, corresponding
            % to no boundary conditions, because it makes no sense imposing a Neumann condition in a vertex»

            if(vertici(i) ~= 1)             % Se il primo vertice non è già di Dirichlet
                vertici(i) = 0;             
            end                             % Altrimenti si lascia invariato
            if(vertici(mod(i,l)+1) ~= 1)    % Se il secondo vertice non è già di Dirichlet
                vertici(mod(i,l)+1) = 0;
            end                             % Altrimenti si lascia invariato
            
        end
    end
    BC.InputVertexValues = vertici;

    
    % Questi indici, oltre a essere marcatori, sono indici che si riferiscono ai valori contenuti nel vettore BC.Values
    % Dunque 1 fa riferimento al primo valore di Dirichlet di BC.Values mentre 6 al sesto di Neumann e cosí via
    
    BC.Holes.Hole = [];
    BC.Segments.Segment = [];
    
    
    
    % --------------------------------------------
    % Inserimento dei parametri di triangolazione
    % --------------------------------------------
    
    RefiningOptions.CheckArea  = 'Y';
    RefiningOptions.CheckAngle = 'N';
    RefiningOptions.AreaValue  = areaMax; % [vale di base: 0.02] Area massima della triangolazione (nessun triangolo avrà area maggiore)
    RefiningOptions.AngleValue = angoloMin; % Minimo angolo consentito in gradianti nei triangoli («must be less then or equal to 30◦ [...]», p. 16 «Manuale BBTR {17-10-2024}»)
    RefiningOptions.Subregions = []; % Finezza controllata a regioni: per raffinare solo alcune parti del dominio
    
    
    % --------------------------------------------
    % Creazione della triangolazione e plottaggio
    % --------------------------------------------
    
    [geom] = mybbtr30(Domain,BC,RefiningOptions);
    % draw_grid (geom,1);
    
    % --------------------------------------------------
    % --------------------------------------------------
    
    
    % --------------------------------------------------
    % Rielaborazione dei prodotti del triangolatore
    % per un piu` agevole trattamento delle condizioni
    % al contorno (operazioni di pulizia)
    % --------------------------------------------------
    
    % Elimina le righe vuote nate dalla generazione
    % (in particolare durante la stima dei triangoli fatta
    % raddoppiando di volta in volta il numero)
    geom.obj.P = geom.obj.P(1:geom.Nobj.N_node,:);
    geom.obj.T = geom.obj.T(1:geom.Nobj.N_ele,:);
    geom.obj.E = geom.obj.E(1:geom.Nobj.N_edge,:);
    geom.obj.Neigh = geom.obj.Neigh(1:geom.Nobj.N_ele,:);
    
    % --------------------------------------------------
    
    % È un micro-«DOFHandler»
    j  = 1;
    Dj = 1;
    for i=1:size(geom.piv.nlist)
         if geom.piv.nlist(i)==0
            geom.piv.piv(i)=j;
            j = j+1;
         else
            geom.piv.piv(i)=-Dj;
            Dj = Dj + 1;
         end
    end
    
    % --------------------------------------------------
    
    geom.piv.piv = transpose(geom.piv.piv);
    
    % --------------------------------------------------
    
    % geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
    % di Dirichlet e il corrispondente marcatore
    
    [X,I] = sort(geom.piv.Di(:,1));
    geom.piv.Di = geom.piv.Di(I,:);
    
    clear X I;

end