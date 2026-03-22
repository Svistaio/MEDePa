function Grafici_Errore(ordP,condA,dimA,errH1,errL2,areeMaxEff,dt,tipoSim)

    % close all

    % In norma H1 si veda eqq. (8.6),(8.7) a p. 92 della filza «Numerical Methods for PDE {23-09-2024}»
    % In norma L2 si veda eqq. (8.16),(8.17) a p. 94 della filza «Numerical Methods for PDE {23-09-2024}»
    % Nota: La scala logaritmica non modifica i valori ma solo le «distanze» tra i valori tramite [ovviamente] il logaritmo

    if(tipoSim == 1) % Se le simulazioni sono in spazio
        xE = sqrt(areeMaxEff);
        tEtiX = "h";
    else
        xE = dt;
        tEtiX = "$\Delta t$";
    end

    if(length(xE)>1)

        xF = linspace(xE(1),xE(end),100);

        if(tipoSim == 1) % Se la simulazione è in spazio

                %%% Numero di condizionamento - Scala logaritmica x-y %%%
            Mostra_Cond(xE,condA,xF)

                %%% Gradi di Libertà - Scala logaritmica x-y %%%
            Mostra_Gdl(xE,dimA,xF)

        end

            %%% Erorri in norma H1 ed L2%%%
        Mostra_Errori(xE,errH1,errL2,xF,ordP,tipoSim,tEtiX)

    else
        disp(['Condizionamento matrice A: ',num2str(condA)])
        disp(['Numero di gdl: ',num2str(dimA)])
        disp(['Errore della soluzione in norma L2: ',num2str(errL2)])
        disp(['Errore della soluzione in norma H1: ',num2str(errH1)])
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funzioni ausiliari %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mostra_Errori(xE,errH1,errL2,xF,ordP,tipoSim,tEtiX)

        global sBottE f tErrLin tErrLog

        sP = []; % Struttura per memorizzare tutte le caratteristiche della figura
        
        f = figure('Name','Errore','NumberTitle','off','Units','normalized', ...
                   'Position',[0 0 1 0.75]); % [left bottom width height]
            movegui(f,'center');

            sSal.stacco = 0.005;
            sSal.dimT = 17.5;
            
                % [left bottom width height]
            sSal.x = 0.2563; sSal.y = 0.05; sSal.l = 0.05; sSal.a = 0.045;
            sSal.P = [sSal.x sSal.y sSal.l sSal.a];
            uicontrol('Style','pushbutton','Enable','on','FontSize',sSal.dimT,'String','Salva', ...
                      'Units', 'Normalized','Position',sSal.P,'CallBack',@FRSalva);

                % [left bottom width height]
            sSal.x = sSal.x+sSal.l+sSal.stacco; sSal.l = 0.1;
            sSal.P = [sSal.x sSal.y sSal.l sSal.a];
   scelta = uicontrol('Style','popupmenu','String',["prima immagine" "seconda immagine" "figura"], ...
                      'Enable','on','FontSize',sSal.dimT,'Units', 'Normalized','Position',sSal.P);
            sBottE.S = scelta;

                % [left bottom width height]
            sSal.x = sSal.x+sSal.l+sSal.stacco;
            sSal.P = [sSal.x sSal.y*7/8 sSal.l sSal.a];
     cntr = uicontrol('Style','text','Enable','on','FontSize',sSal.dimT,'String','in', ...
                      'Units', 'Normalized','Position',sSal.P,'HorizontalAlignment','center');
            e = get(cntr, 'Extent'); sSal.l = e(3);
            cntr.Position(3) = sSal.l;
            % cntr.Position(4)= e(4);

                % [left bottom width height]
            sSal.x = sSal.x+sSal.l+sSal.stacco; sSal.l = 0.17;
            sSal.P = [sSal.x sSal.y sSal.l sSal.a];
            % sSal.nomeImmagine = "";
  percorso = uicontrol('Style','edit','Enable','on','FontSize',sSal.dimT,'String',"./Errore.png", ...
                      'Units', 'Normalized','Position',sSal.P);
            sBottE.P = percorso;

                % [left bottom width height]
            sSal.x = sSal.x+sSal.l+sSal.stacco;
            sSal.P = [sSal.x sSal.y*7/8 sSal.l sSal.a];
     cntr = uicontrol('Style','text','Enable','on','FontSize',sSal.dimT,'String','con risoluzione', ...
                      'Units', 'Normalized','Position',sSal.P,'HorizontalAlignment','center');
            e = get(cntr, 'Extent'); sSal.l = e(3);
            cntr.Position(3) = sSal.l;

                % [left bottom width height]
            sSal.x = sSal.x+sSal.l+sSal.stacco; sSal.l = 0.05;
            sSal.P = [sSal.x sSal.y sSal.l sSal.a];
            sSal.nomeImmagine = "./Triangolazione_Migliore.png";
    risol = uicontrol('Style','edit','Enable','on','FontSize',sSal.dimT,'String','250', ...
                      'Units', 'Normalized','Position',sSal.P);
            sBottE.R = risol;

        dimTria = [0.9,0.7];         % «[lunghezza largezza]»
        hRegione = (1-dimTria)/2; % [altezza regione inferiore, altezza regione superiore]
        tiledlayout(1,2,'TileSpacing','Compact','Padding','loose','Units', 'normalized', ...
                    'Position',[hRegione(1) hRegione(2)+0.05 dimTria(1)+0.03 dimTria(2)]); 

        coloreH1L = 'r';
        coloreH1M = 'k';
        coloreL2L = 'b';
        coloreL2M = 'k';

                %%% Scala lineare x-y %%%
        tErrLin = nexttile;
            
            if(tipoSim == 1)     % Se la simulazione è in spazio
                pH1       = polyfit(xE,errH1,ordP);
                fittLinH1 = polyval(pH1,xF);
                pL2       = polyfit(xE,errL2,(ordP+1));
                fittLinL2 = polyval(pL2,xF);
            elseif(tipoSim == 2) % Se la simulazione è in tempo
                pH1       = polyfit(xE,errH1,2);
                fittLinH1 = polyval(pH1,xF);
                pL2       = polyfit(xE,errL2,2);
                fittLinL2 = polyval(pL2,xF);
            end

            sP.SpessL = 1.5;
            sP.SpessM = 7.5;
            sP.TitT = append("Errore in norma H1 ed L2 dei P",string(ordP)," [scala lin.]");
            sP.TitD = 25;
            sP.EtiD = 17.5;
            sP.LegD = 15;
            cf = 4; % Numero di cifre significative

            Attiva_Griglia(sP.EtiD-5)
            title(tErrLin,sP.TitT,'FontSize',sP.TitD,'Interpreter','latex');

            colororder({coloreH1L,coloreL2L})

            yyaxis left
            yEH1 = errH1; yFH1 = fittLinH1;
                plot(xE,yEH1,'x','Color',coloreH1M,'LineWidth',sP.SpessL,'DisplayName',"Punti dell'errore in $H^1$",'MarkerSize',sP.SpessM);
                hold on
                plot(xF,yFH1,'-','Color',coloreH1L,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio dell'errore in $H^1$");
            alfaLinH1 = pH1(1);KLinH1 = pH1(end);
            t = text(xE(floor(length(xE)/2)+1),yEH1(floor(length(yEH1)/2)+1), ...
                [append("$\alpha^{H1}_{Lin}=$", string(round(alfaLinH1,cf))), ...
                append("$K^{H1}_{Lin}$=",       string(round(KLinH1,cf)))], ...
                'FontSize', 17.5,'Interpreter','latex','Clipping','on',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
            Controllo_Sbordamento(t)

            xlabel(tErrLin,tEtiX,'FontSize',sP.EtiD,'interpreter','latex');
            ylabel(tErrLin,'Errore','FontSize',sP.EtiD,'interpreter','latex','Color','k');
            
            yyaxis right
            yEL2 = errL2; yFL2 = fittLinL2;
                plot(xE,yEL2,'+','Color',coloreL2M,'LineWidth',sP.SpessL,'DisplayName',"Punti dell'errore in $L^2$",'MarkerSize',sP.SpessM);
                plot(xF,yFL2,'--','Color',coloreL2L,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio dell'errore in $L^2$");

            alfaLinL2 = pL2(1); KLinL2 = pL2(end);
            t = text(xE(floor(length(xE)/2)),yEL2(floor(length(yEL2)/2)), ...
                [append("$\alpha^{L2}_{Lin}=$",   string(round(alfaLinL2,cf))), ...
                 append("$K^{L2}_{Lin}$=",        string(round(KLinL2,cf)))], ...
                'FontSize', 17.5,'Interpreter','latex','Clipping', 'on',...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            Controllo_Sbordamento(t)
            ylabel(tErrLin,'','FontSize',sP.EtiD,'interpreter','latex');
            
            axis padded;
            legend(tErrLin, 'Location','northwest','FontSize',sP.LegD,'interpreter','latex');


                %%% scala logaritmica x-y %%%
            tErrLog = nexttile;

            pH1       = polyfit(log(xE),log(errH1),1);
            fittLogH1 = polyval(pH1,log(xF));
            pL2       = polyfit(log(xE),log(errL2),1);
            fittLogL2 = polyval(pL2,log(xF));

            sP.SpessL = 1.5;
            sP.SpessM = 7.5;
            sP.EtiD = 17.5;
            sP.TitT = append("Errore in norma H1 ed L2 dei P",string(ordP)," [scala log.]");
            sP.TitD = 25;
            sP.LegD = 15;
            cf = 4; % Numero di cifre significative


            yEH1 = errH1; yFH1 = exp(fittLogH1);
                loglog(xE,yEH1,'x','Color',coloreH1M,'LineWidth',sP.SpessL,'DisplayName',"Punti dell'errore in $H^1$",'MarkerSize',sP.SpessM);
                hold on
                loglog(xF,yFH1,'-','Color',coloreH1L,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio dell'errore in $H^1$");

            alfaLogH1 = pH1(1); KLogH1 = pH1(2);
            t = text(xE(floor(length(xE)/2)+1),yEH1(floor(length(yEH1)/2)+1), ...
                 [append("$\alpha^{H1}_{Log}=$", string(round(alfaLogH1,cf))), ...
                 append("$K^{H1}_{Log}$=",       string(round(KLogH1,cf)))], ...
                'FontSize', 17.5,'Interpreter','latex','Clipping', 'on',...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            % Controllo_Sbordamento(t)

            yEL2 = errL2; yFL2 = exp(fittLogL2);
                loglog(xE,yEL2,'+','Color',coloreL2M,'LineWidth',sP.SpessL,'DisplayName',"Punti dell'errore in $L^2$",'MarkerSize',sP.SpessM);
                loglog(xF,yFL2,'--','Color',coloreL2L,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio dell'errore in $L^2$");

            alfaLogL2 = pL2(1); KLogL2 = pL2(2);
            t = text(xE(floor(length(xE)/2)+1),yEL2(floor(length(yEL2)/2)+1), ...
                [append("$\alpha^{L2}_{Log}=$",   string(round(alfaLogL2,cf))), ...
                append("$K^{L2}_{Log}$=",         string(round(KLogL2,cf)))], ...
                'FontSize', 17.5,'Interpreter','latex','Clipping', 'on',...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            % Controllo_Sbordamento(t)

            axis padded;
            Attiva_Griglia(sP.EtiD-5)

            title(tErrLog,sP.TitT,'FontSize',sP.TitD,'Interpreter','latex');
            xlabel(tErrLog,tEtiX,'FontSize',sP.EtiD,'interpreter','latex');
            ylabel(tErrLog,'Errore','FontSize',sP.EtiD,'interpreter','latex');
            legend(tErrLog, 'Location','northwest','FontSize',sP.LegD,'interpreter','latex');


        % Scala logaritmica x-y [grafico in scala logaritica ma con valori originali, ove il
        % disegno è corretto e i valori sono leggibili perché non trasformati dal logaritmo stesso]

end

function Mostra_Cond(xE,condA,xF)

        sP = []; % Struttura per memorizzare tutte le caratteristiche della figura

        f = figure('Name','Cond(A)','NumberTitle','off','Units','normalized', ...
                   'Position',[0 0 1 0.75]); % [left bottom width height]
            movegui(f,'center');

        dimTria = [0.9,0.7];         % «[lunghezza largezza]»
        hRegione = (1-dimTria)/2; % [altezza regione inferiore, altezza regione superiore]
        tiledlayout(1,2,'TileSpacing','Compact','Padding','loose','Units', 'normalized', ...
                    'Position',[hRegione(1) hRegione(2)+0.05 dimTria(1)+0.03 dimTria(2)]); 

        coloreF = 'r'; % Colore fittaggio
        coloreP = 'k'; % Colore punti

                %%% Scala lineare x-y %%%
        tErrLin = nexttile;
            
            % pCondA       = polyfit(xE,condA,2);
            % fittLinCondA = polyval(pCondA,xF);

                % Funzione iperbolica
            ramoIp = @(p,h) p(1)./(h.^2)+p(2);
            
                % Errore (somma dei quadrati) della funzione iperbolica
            errRI = @(p) sum((condA-ramoIp(p,xE)).^2);
            
                % Stima iniziale dei parametri
            stimaIP = [1,0];
            
                % Ottimizza i parametri utilizzando «fminsearch»
            pOttimi = fminsearch(errRI,stimaIP);
            
                % Estrai i parametri ottimizzati
            % alfaRI = pOttimi(1);
            % KRI = pOttimi(2);
            
                % Calcola i valori fittati
            fittLinCondA = ramoIp(pOttimi,xF);

            sP.SpessL = 1.5;
            sP.SpessM = 7.5;
            sP.TitT = "Andamento del numero di cond($\mathbf{A}$) [scala lin.]";
            sP.TitD = 25;
            sP.EtiD = 17.5;
            sP.LegD = 15;
            cf = 4; % Numero di cifre significative

            yECondA = condA; %yFCondA = fittLinCondA;
            yFCondA = fittLinCondA;
                plot(xE,yECondA,'x','Color',coloreP,'LineWidth',sP.SpessL,'DisplayName',"Punti di cond($\mathbf{A}$)",'MarkerSize',sP.SpessM);
                hold on
                plot(xF,yFCondA,'-','Color',coloreF,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio di cond($\mathbf{A}$)");
            alfaLinCondA = pOttimi(1);KLinCondA = pOttimi(2);
            t = text(xE(floor(length(xE)/2)+1),yECondA(floor(length(yECondA)/2)+1), ...
                [append("$\alpha^{H1}_{Lin}=$", string(round(alfaLinCondA,cf))), ...
                 append("$K^{H1}_{Lin}$=",       string(round(KLinCondA,cf)))], ...
                 'FontSize', 17.5,'Interpreter','latex','Clipping','on',...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            Controllo_Sbordamento(t)

            axis padded;
            Attiva_Griglia(sP.EtiD-5)

            title(tErrLin,sP.TitT,'FontSize',sP.TitD,'Interpreter','latex');
            xlabel(tErrLin,"h",'FontSize',sP.EtiD,'interpreter','latex');
            ylabel(tErrLin,'cond($\mathbf{A}$)','FontSize',sP.EtiD,'interpreter','latex','Color','k');
            legend(tErrLin, 'Location','northeast','FontSize',sP.LegD,'interpreter','latex');


                %%% scala logaritmica x-y %%%
        tErrLog = nexttile;

            pCondA       = polyfit(log(xE),log(condA),1);
            fittLogCondA = polyval(pCondA,log(xF));

            sP.SpessL = 1.5;
            sP.SpessM = 7.5;
            sP.EtiD = 17.5;
            sP.TitT = "Andamento del numero di cond($\mathbf{A}$) [scala log.]";
            sP.TitD = 25;
            sP.LegD = 15;
            cf = 4; % Numero di cifre significative


            yECondA = condA; yFCondA = exp(fittLogCondA);
                loglog(xE,yECondA,'x','Color',coloreP,'LineWidth',sP.SpessL,'DisplayName',"Punti di cond($\mathbf{A}$)",'MarkerSize',sP.SpessM);
                hold on
                loglog(xF,yFCondA,'-','Color',coloreF,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio di cond($\mathbf{A}$)");

            alfaLogCondA = pCondA(1); KLogCondA = pCondA(2);
            t = text(xE(floor(length(xE)/2)+1),yECondA(floor(length(yECondA)/2)+1), ...
                 [append("$\alpha^{H1}_{Log}=$", string(round(alfaLogCondA,cf))), ...
                  append("$K^{H1}_{Log}$=",       string(round(KLogCondA,cf)))], ...
                  'FontSize', 17.5,'Interpreter','latex','Clipping', 'on',...
                  'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            Controllo_Sbordamento(t)

            axis padded;
            Attiva_Griglia(sP.EtiD-5)

            title(tErrLog,sP.TitT,'FontSize',sP.TitD,'Interpreter','latex');
            xlabel(tErrLog,"h",'FontSize',sP.EtiD,'interpreter','latex');
            ylabel(tErrLog,'cond($\mathbf{A}$)','FontSize',sP.EtiD,'interpreter','latex');
            legend(tErrLog, 'Location','northeast','FontSize',sP.LegD,'interpreter','latex');


        % Scala logaritmica x-y [grafico in scala logaritica ma con valori originali, ove il
        % disegno è corretto e i valori sono leggibili perché non trasformati dal logaritmo stesso]

end

function Mostra_Gdl(xE,dimA,xF)

        sP = []; % Struttura per memorizzare tutte le caratteristiche della figura

        f = figure('Name','Gdl','NumberTitle','off','Units','normalized', ...
                   'Position',[0 0 1 0.75]); % [left bottom width height]
            movegui(f,'center');

        dimTria = [0.9,0.7];         % «[lunghezza largezza]»
        hRegione = (1-dimTria)/2; % [altezza regione inferiore, altezza regione superiore]
        tiledlayout(1,2,'TileSpacing','Compact','Padding','loose','Units', 'normalized', ...
                    'Position',[hRegione(1) hRegione(2)+0.05 dimTria(1)+0.03 dimTria(2)]); 

        coloreF = 'r'; % Colore fittaggio
        coloreP = 'k'; % Colore punti

                %%% Scala lineare x-y %%%
        tErrLin = nexttile;
            
            % pGdl       = polyfit(xE,dimA,2);
            % fittLinGdl = polyval(pGdl,xF);

                % Funzione iperbolica
            ramoIp = @(p,h) p(1)./(h.^2)+p(2);
            
                % Errore (somma dei quadrati) della funzione iperbolica
            errRI = @(p) sum((dimA-ramoIp(p,xE)).^2);
            
                % Stima iniziale dei parametri
            stimaIP = [1,0];
            
                % Ottimizza i parametri utilizzando «fminsearch»
            pOttimi = fminsearch(errRI,stimaIP);
            
                % Estrai i parametri ottimizzati
            % alfaRI = pOttimi(1);
            % KRI = pOttimi(2);
            
                % Calcola i valori fittati
            fittLinGdl = ramoIp(pOttimi, xF);

            sP.SpessL = 1.5;
            sP.SpessM = 7.5;
            sP.TitT = "Andamento del numero di gdl [scala lin.]";
            sP.TitD = 25;
            sP.EtiD = 17.5;
            sP.LegD = 15;
            cf = 4; % Numero di cifre significative

            yEGdl = dimA; % yFGdl = fittLinGdl;
            yFGdl = fittLinGdl;
                plot(xE,yEGdl,'x','Color',coloreP,'LineWidth',sP.SpessL,'DisplayName',"Punti dei gdl",'MarkerSize',sP.SpessM);
                hold on
                plot(xF,yFGdl,'-','Color',coloreF,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio gdl");
            alfaLinH1 = pOttimi(1);KLinH1 = pOttimi(end);
            t = text(xE(floor(length(xE)/2)+1),yEGdl(floor(length(yEGdl)/2)+1), ...
                [append("$\alpha^{H1}_{Lin}=$", string(round(alfaLinH1,cf))), ...
                 append("$K^{H1}_{Lin}$=",       string(round(KLinH1,cf)))], ...
                 'FontSize', 17.5,'Interpreter','latex','Clipping','on',...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            Controllo_Sbordamento(t)

            axis padded;
            Attiva_Griglia(sP.EtiD-5)

            title(tErrLin,sP.TitT,'FontSize',sP.TitD,'Interpreter','latex');
            xlabel(tErrLin,"h",'FontSize',sP.EtiD,'interpreter','latex');
            ylabel(tErrLin,'gdl','FontSize',sP.EtiD,'interpreter','latex','Color','k');
            legend(tErrLin, 'Location','northeast','FontSize',sP.LegD,'interpreter','latex');


                %%% scala logaritmica x-y %%%
        tErrLog = nexttile;

            pGdl       = polyfit(log(xE),log(dimA),1);
            fittLogGdl = polyval(pGdl,log(xF));

            sP.SpessL = 1.5;
            sP.SpessM = 7.5;
            sP.EtiD = 17.5;
            sP.TitT = "Andamento del numero di gdl [scala log.]";
            sP.TitD = 25;
            sP.LegD = 15;
            cf = 4; % Numero di cifre significative


            yEGdl = dimA; yFGdl = exp(fittLogGdl);
                loglog(xE,yEGdl,'x','Color',coloreP,'LineWidth',sP.SpessL,'DisplayName',"Punti dei gdl",'MarkerSize',sP.SpessM);
                hold on
                loglog(xF,yFGdl,'-','Color',coloreF,'LineWidth',sP.SpessL,'DisplayName',"Fittaggio dei gdl");

            alfaLogGdl = pGdl(1); KLogGdl = pGdl(2);
            t = text(xE(floor(length(xE)/2)+1),yEGdl(floor(length(yEGdl)/2)+1), ...
                 [append("$\alpha^{H1}_{Log}=$", string(round(alfaLogGdl,cf))), ...
                  append("$K^{H1}_{Log}$=",       string(round(KLogGdl,cf)))], ...
                  'FontSize', 17.5,'Interpreter','latex','Clipping', 'on',...
                  'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            Controllo_Sbordamento(t)

            axis padded;
            Attiva_Griglia(sP.EtiD-5)

            title(tErrLog,sP.TitT,'FontSize',sP.TitD,'Interpreter','latex');
            xlabel(tErrLog,"h",'FontSize',sP.EtiD,'interpreter','latex');
            ylabel(tErrLog,'gdl','FontSize',sP.EtiD,'interpreter','latex');
            legend(tErrLog, 'Location','northeast','FontSize',sP.LegD,'interpreter','latex');


        % Scala logaritmica x-y [grafico in scala logaritica ma con valori originali, ove il
        % disegno è corretto e i valori sono leggibili perché non trasformati dal logaritmo stesso]

end

function Attiva_Griglia(d)

    grid on % Attiva la griglia nella piastrella
    ax = gca;
    ax.GridLineStyle = '--'; % Cambia le linee della graglia da continue in tratteggiate
    ax.FontSize = d; % Aumenta la dimensione dei caratteri degli assi 

end

function Controllo_Sbordamento(t)

        % Limiti correnti delle ordinate
    yEstr = ylim;

        % Coordinate dell'angolo sinistro del blocco di testo,
        % della sua larghezza e della sua altezza
    estensT = get(t,'Extent');   % Dimensione del testo

        % Se il testo è vicino al limite inferiore, allarga relativo il limite Y
    if(estensT(2)< yEstr(1))
        t.VerticalAlignment = "bottom";
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funzioni di richiamo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = FRSalva(~,~)

    global sBottE f tErrLin tErrLog

    percorso    = string(sBottE.P.String);
    risoluzione = str2double(sBottE.R.String);

    if(sBottE.S.Value == 1)     % Se si vuole salvare la prima immagine
        exportgraphics(tErrLin,percorso,'Resolution',risoluzione)
    elseif(sBottE.S.Value == 2) % Se si vuole salvare la seconda immagine
        exportgraphics(tErrLog,percorso,'Resolution',risoluzione)
    elseif(sBottE.S.Value == 3) % Se si vuole salvare ambo le immagini
        exportgraphics(f,percorso,'Resolution',risoluzione)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%
%%% Codice scartato %%%
%%%%%%%%%%%%%%%%%%%%%%%

    %%% Scala lineare log(x)-log(y) %%%
% p = polyfit(log(h),log(errH1),1);
% fittLog = polyval(p,log(h));
% alfaLog = p(1);
% KLog = p(2);

    % Scala lineare log(x)-log(y)[grafico in scala lineare con valori
    % logaritmici, ove il grafico è corretto ma i valori sono poco
    % leggibili perché trasformati dal logaritmo stesso]
% titoloFig = "Errore in scala lineare log(x)-log(y)";
% titoloGraf = strcat("Errore in norma H1 [scala lineare] - $\alpha_{Log}=$", ...
%                     num2str(alfaLog),", $K_{Log}$=",num2str(KLog));
% x12 = log(h); y1 = log(err); y2 = fittLog;
% creaFigura([x12 x12],[y1 y2],[opzioniGraf1;opzioniGraf2],[titoloFig titoloGraf],["log(h)" "log(ErrH1)"],25,"lin");
