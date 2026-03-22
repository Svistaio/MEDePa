function [] = Visualizza_Triangolazione(ordP,verifica,nN,u_vera,nt,dt,istante,evoluzione)

    % close all

    global geom u_est matrTria
    global zTria c tTria fTria tPunt sBottT
    global p pA pV pD

        % Bottoni di avanzamento, indietreggiamento e selezione istanti degl'istanti temporali
    sBottT = [];

        % Impostazione figura
    f = figure('Name','Triangolazione','NumberTitle','off','Units','Normalized');

        %%% Punti %%%
    if(evoluzione == 1)

        f.Position = [0 0 0.95 0.75];
        % [left bottom width height] e la misura di base sono i fotèli
        % («pixel») normalizzati alla dimensione dello schermo
        movegui(f,'west'); % Muove la figura al sinistra
    
        dimTria = [0.925,0.7];     % «[lunghezza largezza]»
        hRegione = (1-dimTria)/2; % [altezza regione inferiore, altezza regione superiore]

            % Impostazione piastrelle
        tiledlayout(1,2,'TileSpacing','compact', ...
                        'Padding','compact', ...
                        'Units','normalized', ...
                        'Position',[hRegione(1) hRegione(2)+0.05 dimTria(1) dimTria(2)]);
                        % [left bottom width height] e la misura di base sono i fotèli
                        % («pixel») normalizzati alle dimensioni della figura in cui la
                        % piastrella è immersa («f» in questo caso)

        % t.Padding = 'none';
        % t.TileSpacing = 'none';

            % Impostazione triangolazione
        fLinea = nexttile; % Figura in posizione (1,2)

        p = 90;
        X  = (0:nt)*dt;

            yPA = u_est(p,:);
            pA = plot(X,yPA,'-r',"LineWidth",2,'DisplayName','$u_{est}$');

            if(verifica == 1)
                hold on
                
                yPV = u_vera(geom.obj.P(p,1),geom.obj.P(p,2),X);
                pV = plot(X,yPV,'-b',"LineWidth",2,'DisplayName','$u_{vera}$');
        
                yD  = yPV-yPA;
                pD = plot(X,yD,'-k',"LineWidth",2,'DisplayName','diff.');
                hold off
            end

        title('Andamento nel tempo in un nodo/Triangolazione','FontSize',20, 'Units', 'normalized', 'Position', [0.5,1.05],'interpreter','latex');
        xlabel('Tempo','FontSize',15,'interpreter','latex');
        ylabel('u','FontSize',15,'interpreter','latex');
        
        grid on
        ax1 = gca;
        ax1.GridLineStyle = '--'; % Cambia le linee della graglia da continue in tratteggiate

        legend('Location','northeast',... % Posizione della legenda
               'FontSize',17.5,... % Impone la dimensione del carattere nella legenda
               'interpreter','latex'); % Setta l'interpretatore come LaTeX: si possono scrivere formule con quella formattazione all'interno di «$$»


                %%% Bottoni per la navigazione dell'evoluzione in un punto %%%

            % Dimensione bottoni «[left bottom width height]»
        sB.y = 0.05; sB.l = 0.1; sB.a = 0.07;
        sB.stacco = 0.01; dFig = fLinea.Position;
        sB.x = dFig(1)+dFig(3)/2-(3*sB.l+2*sB.stacco)/2-0.005;
        sB.P = [sB.x sB.y sB.l sB.a];
        sB.dimT = 12; % Dimensione testo
        
            % Indientro
       bottIP = uicontrol('Style','pushbutton','Enable','on','FontSize',sB.dimT,'String','Indientro', ...
                          'Units', 'Normalized','Position',sB.P,'CallBack',{@FRIndietroP,u_vera,X,verifica});
        sBottT.IP = bottIP;

            % Nodi
        sB.x = sB.x + sB.l + sB.stacco;
        sB.P = [sB.x sB.y sB.l sB.a];
        listaNodi = append("Nodo ",string(1:nN));
        bottN = uicontrol('Style','popupmenu','Value',p,'Enable','on','FontSize',sB.dimT,'String',listaNodi, ...
                          'Units', 'Normalized','Position',sB.P,'CallBack',{@FRNodi,nN,u_vera,X,verifica});
        sBottT.N = bottN;
        
            % Salto
        sB.a = sB.a/2; sB.l = sB.l/2; sB.x = sB.x + sB.l;
        sB.P = [sB.x sB.y sB.l sB.a];
        bottSP = uicontrol('Style','edit','Units','Normalized','Enable','on','FontSize',sB.dimT,'String','1','Position',sB.P);
        sBottT.SP = bottSP;
        
        sB.x = sB.x - sB.l; sB.P = [sB.x sB.y sB.l sB.a*0.85];
        uicontrol('Style','text','Units','Normalized','Enable','on','FontSize',sB.dimT,'String','Salto:','Position',sB.P);

            % Avanti
        sB.l = sB.l*2; sB.a = sB.a*2; sB.x = sB.x + sB.l + sB.stacco;
        sB.P = [sB.x sB.y sB.l sB.a];
        bottAP = uicontrol('Style','pushbutton','Enable','on','FontSize',sB.dimT,'String','Avanti', ...
                           'Units', 'Normalized','Position',sB.P,'CallBack',{@FRAvantiP,nN,u_vera,X,verifica});
        sBottT.AP = bottAP;

        % FRIstanti(bottN,[],nt)

    else

        f.Position = [0 0 0.5 0.75];
        % [left bottom width height] e la misura di base sono i fotèli
        % («pixel») normalizzati alla dimensione dello schermo

        dimTria = [0.85,0.85];         % «[lunghezza largezza]»
        hRegione = (1-dimTria(1))/2; % [altezza regione inferiore, altezza regione superiore]
    
    
            % Impostazione piastrelle
        tiledlayout(1,1,'TileSpacing','compact', ...
                        'Padding','compact', ...
                        'Units','normalized', ...
                        'Position',[hRegione(1) 0.05 dimTria(1) dimTria(2)]);
                        % [left bottom width height] e la misura di base sono i fotèli
                        % («pixel») normalizzati alle dimensioni della figura in cui la
                        % piastrella è immersa («f» in questo caso)

            % Dimensione bottoni «[left bottom width height]»
        sB.y = 0.05; sB.l = 0.1; sB.a = 0.07;
        sB.stacco = 0.01;
        sB.dimT = 12; % Dimensione testo

    end

    movegui(f,'center'); % Muove la figura al sinistra

    %%% Triangolazione %%%
    c = nt+1;

    Costruisci_Triangolazione(ordP)

            % Impostazione triangolazione
        fTria  = nexttile; % Figura in posizione (1,2)

    if(evoluzione == 1)
        dFig = fTria.Position;
        sB.x = dFig(1)+dFig(3)/2-(3*sB.l+2*sB.stacco)/2+0.005;
        sB.P = [sB.x sB.y sB.l sB.a]; % Si prende la posizione x affinché i tre bottoni siano centrati
    else
        sB.x = 0.5-(3*sB.l+2*sB.stacco)/2;
        sB.P = [sB.x sB.y sB.l sB.a]; % Si prende la posizione x affinché i tre bottoni siano centrati
    end

    zTria = u_est(:,c);
    tTria = trisurf(matrTria,geom.obj.P(:,1),geom.obj.P(:,2),zTria);
    % view(0,0)
    if(evoluzione == 1)
        hold on
        tPunt = scatter3(geom.obj.P(p,1),geom.obj.P(p,2),u_est(p,c),50,'r', 'filled');
        tTria.ButtonDownFcn = @(src,event) Scelta_Punto(src,event,evoluzione,u_vera,X,nN,verifica);
    end
    % title('Triangolazione','FontSize',20, 'Units', 'normalized', 'Position', [0.05,1.05],'interpreter','latex');
    
    grid on
    ax1 = gca;
    ax1.GridLineStyle = '--'; % Cambia le linee della graglia da continue in tratteggiate


    if(evoluzione == 1)
                %%% Bottoni per la navigazione dell'evoluzione della triangolazione %%%
    
            % Indientro
        bottIT = uicontrol('Style','pushbutton','Enable','on','FontSize',sB.dimT,'String','Indientro', ...
                           'Units', 'Normalized','Position',sB.P,'CallBack',{@FRIndietroT,evoluzione});
        sBottT.IT = bottIT;
    
            % Istanti
        sB.x = sB.x + sB.l + sB.stacco;
        sB.P = [sB.x sB.y sB.l sB.a];
        listaIstanti = append("Passo ",string(0:nt)'," - Tempo ",string((0:nt)*dt+istante)');
        bottT = uicontrol('Style','popupmenu','Value',c,'Enable','on','FontSize',sB.dimT,'String',listaIstanti, ...
                          'Units', 'Normalized','Position',sB.P,'CallBack',{@FRIstanti,nt,evoluzione});
        sBottT.T = bottT;
    
            % Salto
        sB.a = sB.a/2; sB.l = sB.l/2; sB.x = sB.x + sB.l;
        sB.P = [sB.x sB.y sB.l sB.a];
        bottST = uicontrol('Style','edit','Units','Normalized','Enable','on','FontSize',sB.dimT, ...
                            'String',num2str(floor(nt/40)),'Position',sB.P);
        sBottT.ST = bottST;
    
        sB.x = sB.x - sB.l; sB.P = [sB.x sB.y sB.l sB.a*0.85];
        uicontrol('Style','text','Units','Normalized','Enable','on',...
                  'FontSize',sB.dimT,'String','Salto:','Position',sB.P);

        sB.P = [0.95 0.01 0.05 0.05];
        bottG = uicontrol('Style','checkbox','Units','Normalized','Enable','on',...
                          'FontSize',sB.dimT,'String','Griglia','Position',sB.P,....
                          'Value',0,'CallBack',{@FRGriglia,f});
        sBottT.G = bottG;

    
            % Avanti
        sB.a = sB.a*2; sB.l = sB.l*2; sB.x = sB.x + sB.l+sB.stacco;
        sB.P = [sB.x sB.y sB.l sB.a];
        bottAT = uicontrol('Style','pushbutton','Enable','on','FontSize',sB.dimT,'String','Avanti', ...
                          'Units', 'Normalized','Position',sB.P,'CallBack',{@FRAvantiT,nt,evoluzione});
        sBottT.AT = bottAT;
    
        FRIstanti(bottT,[],nt,evoluzione)

    else
        sB.P = [0.01 0.95 0.075 0.05];
        bottG = uicontrol('Style','checkbox','Units','Normalized','Enable','on',...
                          'FontSize',sB.dimT,'String','Griglia','Position',sB.P,....
                          'Value',0,'CallBack',{@FRGriglia,f});
        sBottT.G = bottG;
    end

                %%% Bottoni per il Salvataggio %%%

    sSal.stacco = 0.005;
    sSal.dimT = 12.5;
    sSal.RidT = 0.2;
    sSal.l = 0.065;
    sSal.a = 0.045;
    sSal.y = 0.93; 
    if(evoluzione == 1)
        sSal.x = dFig(1)+dFig(3)/2-0.1898+0.005;
    else
        sSal.x = 0.5-0.1898+0.005;
    end
    sSal.P = [sSal.x sSal.y sSal.l sSal.a];
    
        % [left bottom width height]
    uicontrol('Style','pushbutton','Enable','on','FontSize',sSal.dimT,'String','Salva', ...
              'Units', 'Normalized','Position',sSal.P,'CallBack',{@FRSalva,evoluzione});

        % [left bottom width height]
    sSal.x = sSal.x+sSal.l+sSal.stacco;
    sSal.P = [sSal.x sSal.y-sSal.a*sSal.RidT sSal.l sSal.a];
    cntr   = uicontrol('Style','text','Enable','on','FontSize',sSal.dimT,'String','in', ...
                     'Units', 'Normalized','Position',sSal.P,'HorizontalAlignment','center');
    e = get(cntr, 'Extent'); sSal.l = e(3); cntr.Position(3) = sSal.l;

        % [left bottom width height]
    sSal.x   = sSal.x+sSal.l+sSal.stacco; sSal.l = 0.17; 
    sSal.P   = [sSal.x sSal.y sSal.l sSal.a];
    percorso = uicontrol('Style','edit','Enable','on','FontSize',sSal.dimT,'String',"./Triangolazione_Migliore.png", ...
                         'Units', 'Normalized','Position',sSal.P);
    sBottT.P = percorso;

        % [left bottom width height]
    sSal.x = sSal.x+sSal.l+sSal.stacco;
    sSal.P = [sSal.x sSal.y-sSal.a*sSal.RidT sSal.l sSal.a];
    cntr   = uicontrol('Style','text','Enable','on','FontSize',sSal.dimT,'String','con risoluzione ', ...
              'Units', 'Normalized','Position',sSal.P,'HorizontalAlignment','center');
    e = get(cntr, 'Extent'); sSal.l = e(3); cntr.Position(3) = sSal.l;

        % [left bottom width height]
    sSal.x = sSal.x+sSal.l+sSal.stacco; sSal.l = 0.05;
    sSal.P = [sSal.x sSal.y sSal.l sSal.a];
    sSal.nomeImmagine = "./Triangolazione_Migliore.png";
    risol  = uicontrol('Style','edit','Enable','on','FontSize',sSal.dimT,'String','250', ...
              'Units', 'Normalized','Position',sSal.P);
    sBottT.R = risol;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funzioni di richiamo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Triangolazione %%%
function [] = FRIndietroT(~,~,evoluzione)
    
    global c sBottT
    
    salto = str2double(sBottT.ST.String);
    if(c-salto>1)
        c = c-salto;
        Aggiorna_Triangolazione(evoluzione,"on","on")
    else
        c = 1;
        Aggiorna_Triangolazione(evoluzione,"off","on")
    end
    sBottT.T.Value = c;

end

function [] = FRIstanti(srg,~,nt,evoluzione)
    
    global c

    c = srg.Value;
    if(c==nt+1)
        Aggiorna_Triangolazione(evoluzione,"on","off")
    elseif(c==1)
        Aggiorna_Triangolazione(evoluzione,"off","on")
    else
        Aggiorna_Triangolazione(evoluzione,"on","on")
    end

end

function [] = FRAvantiT(~,~,nt,evoluzione)
    
    global c sBottT

    salto = str2double(sBottT.ST.String);
    if(c+salto<nt+1)
        c = c+salto;
        Aggiorna_Triangolazione(evoluzione,"on","on")
    else
        c = nt+1;
        Aggiorna_Triangolazione(evoluzione,"on","off")
    end
    sBottT.T.Value = c;

end

function [] = FRSalva(~,~,evoluzione)

    global fTria sBottT tPunt
    global c u_est matrTria geom

    if(evoluzione == 1)        
        percorso    = string(sBottT.P.String);
        risoluzione = str2double(sBottT.R.String);

        f = figure('Name','Prova','NumberTitle','off','Units','Normalized','Visible','off');
        f.Position = [0 0 0.5 0.75];
        dimTria = [0.85,0.85];         % «[lunghezza largezza]»
        hRegione = (1-dimTria(1))/2; % [altezza regione inferiore, altezza regione superiore]

    
        tiledlayout(1,1,'TileSpacing','compact','Padding','compact','Units','normalized', ...
                        'Position',[hRegione(1) 0.05 dimTria(1) dimTria(2)]);
        t = nexttile;
    
        zTria = u_est(:,c);
        trisurf(matrTria,geom.obj.P(:,1),geom.obj.P(:,2),zTria);
        view(45,45)
        axis tight

        exportgraphics(t,percorso,'Resolution',risoluzione)
        delete(f);

            % In passato lavoravo direttamente sulla figura originale
        % tPunt.Visible = "off";
        % percorso    = string(sBottT.P.String);
        % risoluzione = str2double(sBottT.R.String);
        % % exportgraphics(fTria,percorso,'Resolution',risoluzione)
        % exportgraphics(t,percorso,'Resolution',risoluzione)
        % tPunt.Visible = "on";

    else
        percorso    = string(sBottT.P.String);
        risoluzione = str2double(sBottT.R.String);
        exportgraphics(fTria,percorso,'Resolution',risoluzione)
    end

end

function [] = Aggiorna_Triangolazione(evoluzione,statoI,statoA)

    global u_est
    global zTria c tTria tPunt sBottT
    global p

    zTria = u_est(:,c);
    tTria.Vertices(:,3)   = zTria;    % Aggiorna i punti
    tTria.FaceVertexCData = zTria;    % Aggiorna i colori
    tTria.CDataMapping    = 'scaled'; % Mappatura correta
    drawnow;                          % Aggiorna il grafico

    if(evoluzione == 1)
        tPunt.ZData = u_est(p,c);
    end

    sBottT.IT.Enable = statoI;
    sBottT.AT.Enable = statoA;

    axis tight

end

function [] = Scelta_Punto(~,event,evoluzione,u_vera,X,nN,verifica)

    global geom p sBottT

    if(evoluzione == 1)

            % Coordinate del punto 3D
        pt = event.IntersectionPoint;
        
            % Trova il punto più vicino nella superficie
        [~,p] = min(sqrt((geom.obj.P(:,1)-pt(1)).^2+(geom.obj.P(:,2)-pt(2)).^2));

            % Aggiornamento del punto
        Aggiorna_Punto(u_vera,X,verifica)

        if(p==nN)
            sBottT.AP.Enable = "off";
            sBottT.IP.Enable = "on";
        elseif(p==1)
            sBottT.AP.Enable = "on";
            sBottT.IP.Enable = "off";
        else
            sBottT.AP.Enable = "on";
            sBottT.IP.Enable = "on";
        end

    end
end


    %%% Punti %%%
function [] = FRIndietroP(~,~,u_vera,X,verifica)
    
    global sBottT p
    
    salto = str2double(sBottT.SP.String);

    if(p-salto>1)
        p = p-salto;
        Aggiorna_Punto(u_vera,X,verifica)
        sBottT.AP.Enable = "on";
    else
        p = 1;
        Aggiorna_Punto(u_vera,X,verifica)
        sBottT.IP.Enable = "off";
    end

end

function [] = FRNodi(srg,~,nN,u_vera,X,verifica)
    
    global sBottT p

        p = srg.Value;
        Aggiorna_Punto(u_vera,X,verifica)
        
    if(p==nN)
        sBottT.AP.Enable = "off";
        sBottT.IP.Enable = "on";
    elseif(p==1)
        sBottT.AP.Enable = "on";
        sBottT.IP.Enable = "off";
    else
        sBottT.AP.Enable = "on";
        sBottT.IP.Enable = "on";
    end

end

function [] = FRAvantiP(~,~,nN,u_vera,X,verifica)
    
    global sBottT p

    salto = str2double(sBottT.SP.String);

    if(p+salto<nN)
        p = p+salto;
        Aggiorna_Punto(u_vera,X,verifica)
        sBottT.IP.Enable = "on";
    else
        p = nN;
        Aggiorna_Punto(u_vera,X,verifica)
        sBottT.AP.Enable = "off";
    end

end

function [] = Aggiorna_Punto(u_vera,X,verifica)

        global u_est geom
        global p pA pV pD
        global tPunt c sBottT

        yPA = u_est(p,:);
        pA.YData = yPA;

        tPunt.XData = geom.obj.P(p,1);
        tPunt.YData = geom.obj.P(p,2);
        tPunt.ZData = u_est(p,c);

        if(verifica == 1)
            yPV = u_vera(geom.obj.P(p,1),geom.obj.P(p,2),X);
            pV.YData = yPV;
    
            yD  = yPV-yPA;
            pD.YData = yD;
        end

        sBottT.N.Value = p;

end


    %%% Griglia %%%
function [] = FRGriglia(~,~,f)

    global sBottT geom

    if(isempty(findall(gca,'Type','line')))
        draw_grid(geom,f);
        lines = findall(gca, 'Type', 'line');
        % set(lines,'HitTest','off'); % Disattiva l'interazione col cursore per tutti gli oggetti di tipo «line», ossia la griglia
        % set(lines,'Visible','off'); % Disattiva l'interazione col cursore per tutti gli oggetti di tipo «line», ossia la griglia
    else
        lines = findall(gca, 'Type', 'line');
    end

    if(sBottT.G.Value == 1) % Se il la spunta è premuta s'imposta HitTest e Visible su 'on' per tutte le linee
        set(lines,'Visible','on');
    else % Altrimenti Imposta HitTest e Visible su 'off' per tutte le linee
        set(lines,'Visible','off');
    end

end

%%% %%% %%% %%% %%% %%% %%% 
%%% Funzioni ausiliari  %%% 
%%% %%% %%% %%% %%% %%% %%% 

function Costruisci_Triangolazione(ordP)

    global geom matrTria

        if(ordP == 1)
            matrTria = geom.obj.T;
        elseif(ordP == 2)
            % L'estensione dei P2 fatta dal Berrone estende il vettore degl'indici globali aggiungendo i tre punti intermedi;
            % in particolare il 4° corrisponde a nodo intermedio tra il 1° e il 2°, il 5° tra il 2° e il 3°, infine il 6° tra il 3° e il 1°.
            % Questo giustifica le quattro triadi successive per definire i tre sottotriangoli all'interno di un elemento P2
    
            matrTria = zeros(size(geom.obj.T,1)*4,3);
            for t = 1:size(geom.obj.T,1)
                matrTria((t-1)*4+1,:) = [geom.obj.T(t,1) geom.obj.T(t,4) geom.obj.T(t,6)];
                matrTria((t-1)*4+2,:) = [geom.obj.T(t,2) geom.obj.T(t,5) geom.obj.T(t,4)];
                matrTria((t-1)*4+3,:) = [geom.obj.T(t,3) geom.obj.T(t,6) geom.obj.T(t,5)];
                matrTria((t-1)*4+4,:) = [geom.obj.T(t,6) geom.obj.T(t,4) geom.obj.T(t,5)];
            end % for t
        else
            disp("L'ordine inserito non è contemplato.")
        end

            % Esporta la figura «f» con nome «Nome_Filza» e formato «formato» (a.e. «.png», «.jpg», ecc.) nel percorso «Percoso» (a.e. «C:/Users/Valerio/Desktop/Scaricati»)
        % exportgraphics(f,'C:/Users/Valerio/Desktop/Scaricati/DNO P2 0.005.png',)

end


