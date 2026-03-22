function Evolutore(evoluzione,u0,dt,nt,istante,ngdlMax,gD)

    global geom A D C M b Ad Dd Cd Md ud u u_est

    if(evoluzione == 0) % Se il problema non è evolutivo 

            % Calcolo del vettore di Dirichlet
        ud(:) = gD(geom.obj.P(geom.piv.Di(:,1),1),geom.obj.P(geom.piv.Di(:,1),2),istante);

            % Calcolo della soluzione
        u = A\(b-Ad*ud);
    
            % Estensione della soluzione per considerare le condizioni di Dirichlet
        u_est = zeros(geom.Nobj.N_node,nt+1); % Soluzione estesa ai punti di Dirichlet
        Estendi_Soluzione(1)

    else % Altrimenti

            % Inizializzazioni
        u   = zeros(ngdlMax,nt+1);
        ud  = zeros(size(ud,1),nt+1);  % Vettore di Dirichlet
            % Si considera nt+1 perché la prima colonna è la soluzione iniziale

            % Condizioni iniziale
        for i = 1:length(geom.piv.piv)
            ii = geom.piv.piv(i);
            if(ii>0)
                u(ii,1) = u0(geom.obj.P(i,1),geom.obj.P(i,2));
            else % Condizione iniziale del vettore di Dirichlet
                ud(-ii,1)  = gD(geom.obj.P(geom.piv.Di(-ii,1),1),geom.obj.P(geom.piv.Di(-ii,1),2),0*dt+istante);
            end
        end
        
            % Estensione della soluzione all'istante iniziale
        u_est = zeros(geom.Nobj.N_node,nt+1); % Soluzione estesa ai punti di Dirichlet
        Estendi_Soluzione(1)

            % Calcolo della matrice evolutiva
        Ae = M+(dt/2)*A; % Matrice evolutiva
        
        for n=1:nt
                % Comunicazione utente
            disp(append("Calcolo del ",string(n),"° passo."))

                % Calcolo dei vettori temporali al tempo successivo
            ud(:,n+1)  = gD(geom.obj.P(geom.piv.Di(:,1),1),geom.obj.P(geom.piv.Di(:,1),2),n*dt+istante);

                % Ricalcolo del termine noto evolutivo
            be = (M-(dt/2)*A)*u(:,n)-(dt/2)...
                  *((2/dt)*Md*(ud(:,n+1)-ud(:,n)) ...
                  + Ad*(ud(:,n+1)+ud(:,n)) ...
                  - (b(:,n+1)+b(:,n))); % Termine noto evolutivo

                % Soluzione al passo (n+1)-esimo
            u(:,n+1) = Ae\be;

                % Estensione delle soluzione in tutti gl'istanti di tempo per considerare le condizioni di Dirichlet
            Estendi_Soluzione(n+1)
        end

        % Il termine «(2/dt)*Md*(ud(:,n+1)-ud(:,n))» è estremamente
        % peculiare e mi ha tormentato non poco: in pratica, rispetto alla
        % teoria vista a lezione (p. 51 Appunti Berrone 2024), approssima
        % all'istante t^n la derivata u'_d(t^n) in avanti mentre quella
        % u'_d(t^(n+1)) all'indietro

        % In altre parole in tal modo u'_d assume due valori della derivata
        % in t^n: all'instante t^n quello della derivata in avanti, mentre
        % all'istante t^(n+1) quello della derivata all'indietro, le cui
        % approssimazioni possono potenzialmente differire fra di loro

        % Ora il perché si faccia questo il professore Berrone non lo ha
        % minimamente spiegato, limitandosi a lezione di spiegare che le
        % derivate si approssimavano mediante differenze finite.
        
        % Se dovessi dare una motivazione, è probabile che in tal modo si
        % evita la necessità di avere due condizioni iniziali perché,
        % infatti, gli elementi sono sempre due: n ed n+1. Tuttavia ricordo
        % molto chiaramente che, provando ad approssimare ambo le derivate
        % colla medesima differenza finita, non riuscivo a ottenere
        % un risultato perfetto, essendo presente sempre un lieve errore.

        % Eppure credo d'aver trovato la ragione: è per coerenza dei metodi
        % usati; infatti il metodo di Crank-Nicolson non è che una media tra
        % il metodo di Eulero esplicito e implicto, e questi due approssimano
        % u'_d(t^n) viene dalla parte esplicita mentre u'_d(t^(n+1)) da quella
        %implicita; però queste, per valere, approssimano le derivate temporali
        % in avanti e all'indietro rispetto al tempo valutato [1, 2]. Ecco perché
        % approssimando u'_d(t^n) all'indietro e u'_d(t^(n+1)) in avanti si
        % commetteva un lieve errore dovuto alla semplice incoerenza coi metodi di partenza

        % [1] Dd. 390 e 394 della filza «Metodi Numerici {27-09-2022}.pdf»
        % [2] https://en.wikipedia.org/wiki/Explicit_and_implicit_methods#Illustration_using_the_forward_and_backward_Euler_methods

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funzioni ausiliari %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Estendi_Soluzione(n)

    global geom ud u u_est

    for i=1:geom.Nobj.N_node
        ii = geom.piv.piv(i);
        if(ii>0)
            u_est(i,n) = u(ii,n);
        else
            % u_est(i) = gD(geom.obj.P(geom.piv.Di(-ii),1),geom.obj.P(geom.piv.Di(-ii),2));
            u_est(i,n) = ud(-ii,n);
        end
    end

end