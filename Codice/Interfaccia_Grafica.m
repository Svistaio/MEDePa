function [datiS,datiE,datiD,datiF,datiT] = Interfaccia_Grafica()

    % Variabili globali per tracciare se sono stati premuti i bottoni; infatti Matlab, per qualche motivazione,
    % cambia «AvvioBott.Value» e «EsciBott.Value» da 0 a 1 solo finché si è all'interno della funzione di richiamo:
    % una volta fuori il suo valore è riportato a zero, come se non fosse stato premuto
    % (si v. «https://it.mathworks.com/matlabcentral/answers/394250-how-do-i-define-a-pushbutton-to-change-the-value-of-a-variable»)
    global avvia termina
    avvia = 0; termina = 0;

%%
    % Matrici contenenti in modo estremamente sintetico le funzioni e i casi di studio rispettivamente
    global sDom sParam sFunz sCS


        % Casi di studio
    % In totale vi sono 29 parametri da fornire al programma (ove alcune informazioni sono relative alla geometria)
    sCS = []; sCS.NCS = "Caso di Studio all'avvio";

    syms x y
    v = [x y];

        n = "CSA"; sCS.(n).Nome = "Caso di Studio all'avvio";
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 1;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "10^-6"; sCS.(n).Beta = "[1 0]";   sCS.(n).Sigma = "0";
    sCS.(n).Area_Max = "0.01"; sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    f = 16*x*(1-x)*y*(1-y)*exp(4*x*y)+x+y;
    sCS.(n).Rifu  = strrep(string(f),"*",".*");
    gf = gradient(f,v); gf = string(gf);
    gf = strrep(gf,"*",".*"); gf = strrep(gf,"/","./"); gf = strrep(gf,"^",".^");
    sCS.(n).Rifu  = strrep(string(f),"*",".*");
    sCS.(n).dRifu = [gf(1),gf(2),"0"];
    lf = laplacian(f,v); lf = string(lf);
    lf = strrep(lf,"*",".*"); lf = strrep(lf,"/","./"); lf = strrep(lf,"^",".^");
    sCS.(n).Riff  = string(lf);
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        % Funzioni interessanti da studiare
    % A=10;
    % B=15;
    % sCS.(n).Rifu  = append(string(A),"*exp(",string(-B),"*x).*exp(",string(-B),"*y)");
    % sCS.(n).dRifu = [append(string(-A*B),"*exp(",string(-B),"*x).*exp(",string(-B),"*y)"),...
    %                  append(string(-A*B),"*exp(",string(-B),"*x).*exp(",string(-B),"*y)"),...
    %                  "0"];
    % sCS.(n).Riff  = append(string(A*2*B^2),"*exp(",string(-B),"*x).*exp(",string(-B),"*y)");


    % sCS.(n).dRifu = ["exp(x+y)+2.*(1+x+y)+cos(x).*cos(pi.*x).*cos(pi.*y).*sin(y)-pi.*cos(pi.*y).*sin(x).*sin(pi.*x).*sin(y)", ...
    %                  "exp(x+y)+2.*(1+x+y)+cos(pi.*x).*sin(x).*(cos(y).*cos(pi.*y)-pi.*sin(y).*sin(pi.*y))",...
    %                  "0"];
    % sCS.(n).Riff  = "4+2.*exp(x+y)-2.*cos(pi.*x).*cos(pi.*y).*sin(x).*sin(y)-2.*(pi.^2).*cos(pi.*x).*cos(pi.*y).*sin(x).*sin(y)" + ...
    %                 "-2.*pi.*cos(x).*cos(pi.*y).*sin(pi.*x).*sin(y)-2.*pi.*cos(pi.*x).*cos(y).*sin(x).*sin(pi.*y)";
    
      % "5*sin(2*pi.*y)*cos(2*pi.*x)"
      % "-10.*pi.*sin(2*pi.*x).*sin(2*pi.*y)" "10.*pi.*cos(2*pi.*x).*cos(2*pi.*y)"
      % "-20.*pi.^2.*cos(2*pi.*x).*sin(2*pi.*y)"

      % "1./((x+1).*(y+1))-exp(-((x-0.1).^2+(y-0.1).^2)./(2*0.01))"
      % "100.*exp(-50.*((-0.1+x).^2+(-0.1+y).^2)).*(-0.1+x)-1./((1+x).^2.*(1+y))" "100.*exp(-50.*((-0.1+x).^2+(-0.1+y).^2)).*(-0.1+y)-1/((1+x).*(1+y).^2)"
      % "exp(-50.*(-0.1+x).^2-50.*(-0.1+y).^2).*(-2.84217.*10.^(-14)+2000.*.Intervallox-10000.*x.^2+2000.*y-10000.*y.^2)+(2.*(2+2.*x+x.^2+2.*y+y.^2))./((1+x).^3.*(1+y).^3)"

            % Casi di studio classici
        n = "PM1"; sCS.(n).Nome = "Precisione di macchina P1"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "x+y";
    sCS.(n).dRifu = ["1", "1" "0"];
    sCS.(n).Riff  = "0";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "PM2"; sCS.(n).Nome = "Precisione di macchina P2"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "x.^2+y.^2";
    sCS.(n).dRifu = ["2*x", "2*y" "0"];
    sCS.(n).Riff  = "4";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "DO"; sCS.(n).Nome = "Sola Diffusione DO"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[0 0]";   sCS.(n).Sigma = "0";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "16*x.*(1-x).*y.*(1-y)";
    sCS.(n).dRifu = ["16.*(1-2*x).*y.*(1-y)", "16.*(1-2*y).*x.*(1-x)" "0"];
    sCS.(n).Riff  = "32*(x.*(x-1)+y.*(y-1))";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "DNO"; sCS.(n).Nome = "Sola Diffusione DNO"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 2; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[0 0]";   sCS.(n).Sigma = "0";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "16*x.*(1-x).*y.*(1-y)+x+y";
    sCS.(n).dRifu = ["16.*(1-2*x).*y.*(1-y)+1", "16.*(1-2*y).*x.*(1-x)+1" "0"];
    sCS.(n).Riff  = "32*(x.*(x-1)+y.*(y-1))";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "SDC1"; sCS.(n).Nome = "Stabilità Diffusione-Convezione I"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 1;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "10^-6"; sCS.(n).Beta = "[1 0]";   sCS.(n).Sigma = "0";
    sCS.(n).Area_Max = "0.01"; sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "16*x.*(1-x).*y.*(1-y).*exp(4*x)+x+y";
    sCS.(n).dRifu = ["1+16.*exp(4.*x).*(-1-2.*x+4.*x.^2).*(-1+y).*y", "1+16.*exp(4.*x).*(-1+x).*x.*(-1+2.*y)" "0"];
    sCS.(n).Riff  = "32.*exp(4.*x).*(-x-3.*(-1+y).*y+x.^2.*(1-8.*y+8.*y.^2))";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;
    
        n = "SDC2"; sCS.(n).Nome = "Stabilità Diffusione-Convezione II"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 1;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "10^-6"; sCS.(n).Beta = "[1 0]";   sCS.(n).Sigma = "0";
    sCS.(n).Area_Max = "0.01"; sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    alfa = 0.005; %0.005 Rapidità di decadimento: può variare da 0 [escluso] in poi e piú è piccola piú la funzione raassomiglia a un cubo di ampiezza unitaria
    epsilon = 0.0005; % Costante per garantire la liscezza delle derivate al bordo
    norm = exp(alfa/(0.5^2/(1+epsilon)+epsilon)^2); % Costante di normalizzazione
    f = norm*exp(-alfa/( ...
             ((x*(1-x))/(1+epsilon)+epsilon) ...
             *((y*(1-y))/(1+epsilon)+epsilon)...
             ));
    sCS.(n).Rifu  = strrep(strrep(strrep(string(f),"*",".*"),"/","./"),"^",".^");
    gf = gradient(f,v); gf = string(gf);
    gf = strrep(gf,"*",".*"); gf = strrep(gf,"/","./"); gf = strrep(gf,"^",".^");
    sCS.(n).dRifu = [gf(1),gf(2),"0"];
    lf = laplacian(f,v); lf = string(lf);
    lf = strrep(lf,"*",".*"); lf = strrep(lf,"/","./"); lf = strrep(lf,"^",".^");
    sCS.(n).Riff  = lf;
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;
    
        n = "SDR1"; sCS.(n).Nome = "Stabilità Diffusione-Reazione I"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 1;
    sCS.(n).Ni = "10^-6"; sCS.(n).Beta = "[0 0]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    a = 10;
    f1 = exp(a*(x-1)).*exp(a*(y-0.5))+exp(a*(x-1)).*exp(-a*(y-0.5));
    f2  = exp(-a*x)*exp(-a*(y-0.5))+exp(-a*x)*exp(a*(y-0.5));
    fts = f1*f2; fti = subs(f1,[x,y],[y,x])*subs(f2,[x,y],[y,x]);
    f = -fts-fti+1;
    gf = gradient(f,v); gf = string(gf);
    gf = strrep(gf,"*",".*"); gf = strrep(gf,"/","./"); gf = strrep(gf,"^",".^");
    sCS.(n).Rifu  = strrep(string(f),"*",".*");
    sCS.(n).dRifu = [gf(1),gf(2),"0"];
    lf = laplacian(f,v); lf = string(lf);
    lf = strrep(lf,"*",".*"); lf = strrep(lf,"/","./"); lf = strrep(lf,"^",".^");
    sCS.(n).Riff  = string(lf);
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "SDR2"; sCS.(n).Nome = "Stabilità Diffusione-Reazione II"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 1;
    sCS.(n).Ni = "10^-6"; sCS.(n).Beta = "[0 0]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    alfa = 0.005; % Rapidità di decadimento: può variare da 0 [escluso] in poi e piú è piccola piú la funzione raassomiglia a un cubo di ampiezza unitaria
    epsilon = 0.0005; % Costante per garantire la liscezza delle derivate al bordo
    norm = exp(alfa/(0.5^2/(1+epsilon)+epsilon)^2); % Costante di normalizzazione
    f = norm*exp(-alfa/( ...
             ((x*(1-x))/(1+epsilon)+epsilon) ...
             *((y*(1-y))/(1+epsilon)+epsilon)...
             ));
    sCS.(n).Rifu  = strrep(strrep(strrep(string(f),"*",".*"),"/","./"),"^",".^");
    gf = gradient(f,v); gf = string(gf);
    gf = strrep(gf,"*",".*"); gf = strrep(gf,"/","./"); gf = strrep(gf,"^",".^");
    sCS.(n).dRifu = [gf(1),gf(2),"0"];
    lf = laplacian(f,v); lf = string(lf);
    lf = strrep(lf,"*",".*"); lf = strrep(lf,"/","./"); lf = strrep(lf,"^",".^");
    sCS.(n).Riff  = lf;
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "PME1"; sCS.(n).Nome = "Precisione di macchina P1-CN"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[0 0]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "x+y+t.^2";
    sCS.(n).dRifu = ["1", "1" "2*t"];
    sCS.(n).Riff  = "0";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 1; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 2;

        n = "PME2"; sCS.(n).Nome = "Precisione di macchina P2-CN"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 2; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[0 0]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "x.^2+y.^2+t.^2";
    sCS.(n).dRifu = ["2*x", "2*y" "2*t"];
    sCS.(n).Riff  = "4";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 1; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 2;

        n = "ET1"; sCS.(n).Nome = "Problema Evolutivo - P1-CN (Errore in Tempo)"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "x+y+sin(5*t)";
    sCS.(n).dRifu = ["1", "1" "cos(5*t)*5"];
    sCS.(n).Riff  = "0";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 1; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 2;

        n = "ET2"; sCS.(n).Nome = "Problema Evolutivo - P2-CN (Errore in Tempo)"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 2; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "x.^2+y.^2+sin(5*t)";
    sCS.(n).dRifu = ["2*x", "2*y" "cos(5*t)*5"];
    sCS.(n).Riff  = "4";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 1; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 2;

        n = "ES12"; sCS.(n).Nome = "Problema Evolutivo - P1/P2-CN (Errore in Spazio)"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 1; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "sin(x)+y.^3+t.^2";
    sCS.(n).dRifu = ["cos(x)", "3*y.^2" "2*t"];
    sCS.(n).Riff  = "-sin(x)+6*y";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 1; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "EGM"; sCS.(n).Nome = "Problema Evolutivo - Gaussiana Movente"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 1; sCS.(n).OrdineP = 2; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "3";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 2 2 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "0.25*exp(-20*((x-t).^2+(y-t).^2))";
    sCS.(n).dRifu = ["0.25*exp(-20*((x-t).^2+(y-t).^2))*(-40*(x-t))", ...
                     "0.25*exp(-20*((x-t).^2+(y-t).^2))*(-40*(y-t))", ...
                     "0.25*exp(-20*((x-t).^2+(y-t).^2))*(40*(x+y-2*t))"];
    sCS.(n).Riff  = "0.25*exp(-20*((x-t).^2+(y-t).^2))*(1600*((x-t).^2+(y-t).^2)-80)";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 1; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;

        n = "DOSV"; sCS.(n).Nome = "Sola Diffusione DNO Senza Verifica"; sCS.NCS(end+1) = sCS.(n).Nome;
    sCS.(n).Verifica = 0; sCS.(n).OrdineP = 2; sCS.(n).Stabilita = 0;  sCS.(n).Ammassamento = 0;
    sCS.(n).Ni = "1"; sCS.(n).Beta = "[1 1]";   sCS.(n).Sigma = "1";
    sCS.(n).Area_Max = "0.01";     sCS.(n).Theta_min = "20"; sCS.(n).NSimulazioni = "5";  sCS.(n).Base = "2";
    sCS.(n).Nome_Dominio = "Quadrato"; sCS.(n).Info_Contorno = ["Bordi" "[1 2 2 1]"];
    sCS.(n).Car1 = ["Lato" "on" "1"]; sCS.(n).Car2 = ["Angolo" "on" "[0 0]"]; sCS.(n).Car3 = ["Vuoto" "off" "vuoto"];
    sCS.(n).Rifu  = "16*x.*(1-x).*y.*(1-y)";
    sCS.(n).dRifu = ["16.*(1-2*x).*y.*(1-y)", "16.*(1-2*y).*x.*(1-x)" "16*x.*(1-x).*y.*(1-y)"];
    sCS.(n).Riff  = "-32*(x.*(x-1)+y.*(y-1))+16.*(1-2*x).*y.*(1-y)+16.*(1-2*y).*x.*(1-x)+16*x.*(1-x).*y.*(1-y)";
    sCS.(n).DeltaT  = "0.01"; sCS.(n).Intervallo  = "1"; sCS.(n).NaturaP  = 1;
    sCS.(n).Evoluzione = 0; sCS.(n).Istante = "0"; sCS.(n).TipoSim = 1;


        % Domini
    sDom = []; sDom.NDom = "Nomi dei domini";

        n = "Quadrato"; sDom.(n).Nome = "Dominio a Quadrato"; sDom.NDom(end+1) = sDom.(n).Nome;
    sDom.(n).Valore_Dominio = 1; sDom.(n).Info_Contorno = ["Bordi" "[1 1 1 1]"];
    sDom.(n).Car1 = ["Lato" "on" "1"]; sDom.(n).Car2 = ["Angolo" "on" "[0 0]"]; sDom.(n).Car3 = ["Vuoto" "off" "vuoto"];

        n = "Elle"; sDom.(n).Nome = "Dominio a Elle"; sDom.NDom(end+1) = sDom.(n).Nome;
    sDom.(n).Valore_Dominio = 2; sDom.(n).Info_Contorno = ["Bordi" "[1 1 1 1 1 1]"];
    sDom.(n).Car1 = ["Lato" "on" "1"]; sDom.(n).Car2 = ["Angolo" "on" "[0 0]"]; sDom.(n).Car3 = ["Vuoto" "off" "vuoto"];

        n = "Triangolo"; sDom.(n).Nome = "Dominio a Triangolo"; sDom.NDom(end+1) = sDom.(n).Nome;
    sDom.(n).Valore_Dominio = 3; sDom.(n).Info_Contorno = ["Bordi" "[1 1 1]"];
    sDom.(n).Car1 = ["rad(Area*2)" "on" "1"]; sDom.(n).Car2 = ["Lambda" "on" "0.1"]; sDom.(n).Car3 = ["Base" "on" "[0 0]"];

        n = "Cerchio"; sDom.(n).Nome = "Dominio a Cerchio"; sDom.NDom(end+1) = sDom.(n).Nome;
    sDom.(n).Valore_Dominio = 4; sDom.(n).Info_Contorno = ["% Dirichlet" "1"];
    sDom.(n).Car1 = ["N. punti" "on" "100"]; sDom.(n).Car2 = ["Raggio" "on" "sqrt(2)/2"]; sDom.(n).Car3 = ["Centro" "on" "[0 0]"];


        % Parametri
    sParam = []; sParam.NParam = "Nomi delle configurazioni dei parametri";

        n = "SD"; sParam.(n).Nome = "Sola Diffusione"; sParam.NParam(end+1) = sParam.(n).Nome;
    sParam.(n).Ni = "1"; sParam.(n).Beta = "[0 0]";   sParam.(n).Sigma = "0";
        n = "DC"; sParam.(n).Nome = "Diffusione-Convezione"; sParam.NParam(end+1) = sParam.(n).Nome;
    sParam.(n).Ni = "10^-6"; sParam.(n).Beta = "[1 0]";   sParam.(n).Sigma = "0";
        n = "DR"; sParam.(n).Nome = "Diffusione-reazione"; sParam.NParam(end+1) = sParam.(n).Nome;
    sParam.(n).Ni = "10^-6"; sParam.(n).Beta = "[0 0]";   sParam.(n).Sigma = "1";


        % Funzioni
    sFunz = []; sFunz.NFunz = "Nomi delle funzioni";

        n = "V"; sFunz.(n).Nome = "Vuoto"; sFunz.NFunz(end+1) = sFunz.(n).Nome;
    sFunz.(n).Verifica = 1; sFunz.(n).Rifu = "0"; sFunz.(n).dRifu = ["0" "0"]; sFunz.(n).Riff = "0";

        n = "OU05"; sFunz.(n).Nome = "Omogeneo e unitario in (0.5,0.5)"; sFunz.NFunz(end+1) = sFunz.(n).Nome;
    sFunz.(n).Verifica = 1; sFunz.(n).Rifu = "16*x.*(1-x).*y.*(1-y)";
    sFunz.(n).dRifu = ["16.*(1-2*x).*y.*(1-y)", "16.*(1-2*y).*x.*(1-x)"];
    sFunz.(n).Riff = "32*(x.*(1-x)+y.*(1-y))";

        n = "OU0"; sFunz.(n).Nome = "Omogeneo e unitario all'origine"; sFunz.NFunz(end+1) = sFunz.(n).Nome;
    sFunz.(n).Verifica = 1; sFunz.(n).Rifu = "16.*(1/2+x).*(1/2-x).*(1/2+y).*(1/2-y)";
    sFunz.(n).dRifu = ["16.*(1/2-x).*(1/2+y).*(1/2-y)-16.*(1/2+x).*(1/2+y).*(1/2-y)", "16.*(1/2-y).*(1/2+x).*(1/2-x)-16.*(1/2+y).*(1/2+x).*(1/2-x)"];
    sFunz.(n).Riff = "-16.*(-1+2.*x.^2+2.*y.^2)";
    
        n = "RG"; sFunz.(n).Nome = "Rampa a gobba"; sFunz.NFunz(end+1) = sFunz.(n).Nome;
    sFunz.(n).Verifica = 1; sFunz.(n).Rifu = "16*x.*(1-x).*y.*(1-y)+x+y";
    sFunz.(n).dRifu = ["16.*(1-2*x).*y.*(1-y)+1", "16.*(1-2*y).*x.*(1-x)+1"];
    sFunz.(n).Riff = "32*(x.*(1-x)+y.*(1-y))";

        n = "ACaso"; sFunz.(n).Nome = "Prova di valori un po' a caso: si gioca!"; sFunz.NFunz(end+1) = sFunz.(n).Nome;
    sFunz.(n).Verifica = 0; sFunz.(n).Rifu = "16*x.*(1-x).*y.*(1-y)";
    sFunz.(n).dRifu = ["1", "1"];
    sFunz.(n).Riff = "-x.^3.*y.^3";


%%

        % Posizione e dimensione della figura (bella anche la configurazione [0.1 0.1 0.8 0.8])
    pf = [0 0];     % Posizione della figura
    df = [1/1.9 1/1.4]; % Dimensione della figura 
        % «[left bottom width height]»

    figIn = figure('Name','Interfaccia Grafica', ...
                   'Units','Normalized', ...
                   'Position',[pf(1) pf(2) df(1) df(2)], ... % In origine v'era «[1250 300 560 420]»
                   'Resize','off', ... % Dà («on») o toglie («off») la possibilità di ridimensionare la figura
                   'NumberTitle','off');
        
        movegui(figIn,'center')

        % Funzione «uicontrol» di Matlab, si v. https://it.mathworks.com/help/matlab/ref/uicontrol.html
        % per l'argomento «'Style'» si v. https://it.mathworks.com/help/matlab/ref/uicontrol.html#namevaluepairarguments

            % Spiegazione della funzione e delle relative componenti
        % uicontrol('Style','text',...                                 % Oggetto a testo statico
        %           'Style','edit',...                                 % Oggetto a testo dinamico
        %           'Style','popup',...                                % Oggetto a tendina   
        %           'FontSize',12,...                                  % Dimensione testo
        %           'String',testoCol(1),...                           % Stringa nell'oggetto
        %           'Value',1,...
        %           'HorizontalAlignment','center',...                 % Allineamento del testo
        %           'Units', 'Normalized',...                          % Unità normalizzate
        %           'Position',[xCol(1) yRig(1) dOggParam(1) dOggParam(2)],...   % Posizione «[left bottom width height]»
        %           'Visible','on',...                                 % Visibilità dell'oggetto
        %           'Enable','on');                                    % Abilitazione all'interazione

        global oggIG
        oggIG = []; % Struttura per memorizzare tutti gli oggetti dell'interfaccia grafica


            % Ascisse delle colonne di riferimento
        colRif      = linspace(0.25,0.75,3);  % «[left]»
        staccoOgg   = 0.025;

            % Ordinate della righe di riferimento
        rigRif      = [0.85 linspace(0.75,0.3,6)]; % «[bottom]»
        staccoTesto = 0.03; % Stacco inferiore tra testo statico e dinamico

            % Dimensione delle caselle dinamiche e statiche
        dOggParam = [.1 .04]; % «[width height]»
        dOggFDSel = [.15 .04]; % «[width height]»
        dOggFunz  = [.45 .04]; % «[width height]»

        testoOgg  = [];
        testoOgg.OrdP = "Ordine P"; testoOgg.Stab = "Stabilità";
        testoOgg.Verifica = "Verifica"; testoOgg.Ammassamento = "Massa";
        testoOgg.Ni = "Ni"; testoOgg.Beta = "Beta"; testoOgg.Sigma =  "Sigma";
        testoOgg.SelParam = "Selezione parametri"; testoOgg.SelDom = "Selezione Dominio";
        testoOgg.SelFunz = "Selezione funzioni"; testoOgg.CS = "Casi di studio";
        testoOgg.Area_Max = "A_max"; testoOgg.Theta_min = "θ_min";
        testoOgg.NSimulazioni = "N. Sim."; testoOgg.Base = "Base";
        testoOgg.Dt = "Δt"; testoOgg.Intervallo = "Intervallo T"; testoOgg.NaturaP = "Natura Eq.";
        testoOgg.Evoluzione = "Evoluzione";

        dimTOgg = 10;        % Dimensione testo

%%

            % % % Commenti generali % % %
                    uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','Nota: per il campo «Bordi» vale Dirichlet==1 e Neumann==2','HorizontalAlignment','left',...
                              'Units', 'Normalized','Position',[0.01 0.95 1 dOggParam(2)],...
                              'Visible','on','Enable','on');
                    uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','Nota: non porre Neumann nei lati con β entrante: si perde la coercività e la soluzione esplode','HorizontalAlignment','left',...
                              'Units', 'Normalized','Position',[0.01 0.92 1 dOggParam(2)],...
                              'Visible','on','Enable','on');
                    uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',"Nota: l'unico problema evolutivo considerato è parabolico.",'HorizontalAlignment','left',...
                              'Units', 'Normalized','Position',[0.0875 0.25 1 dOggParam(2)],...
                              'Visible','on','Enable','on');


            % % % Prima colonna % % %
         testoOrd = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.OrdP,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(2)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
          datoOrd = uicontrol('Style','popupmenu','Value',sCS.CSA.OrdineP,...
                              'String',["1" "2"],'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(2) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.OrdineP = [datoOrd testoOrd];
    
          testoNi = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Ni,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(3)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
           datoNi = uicontrol('Style','edit',...
                              'String',sCS.CSA.Ni,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(3) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.Ni = [datoNi testoNi];
    
        testoBeta = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Beta,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(4)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
         datoBeta = uicontrol('Style','edit',...
                              'String',sCS.CSA.Beta,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(4) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.Beta = [datoBeta testoBeta];
    
       testoSigma = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Sigma,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(5)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        datoSigma = uicontrol('Style','edit',...
                              'String',sCS.CSA.Sigma,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(5) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.Sigma = [datoSigma testoSigma];


            % % % Quinta colonna % % %
        testoLati = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',sCS.CSA.Info_Contorno(1),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(2)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
         datoLati = uicontrol('Style','edit',...
                              'String',sCS.CSA.Info_Contorno(2),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(2) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.Lati = [datoLati testoLati];

        testoCar1 = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',sCS.CSA.Car1(1),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(3)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible',sCS.CSA.Car1(2),'Enable','on');
         datoCar1 = uicontrol('Style','edit',...
                              'String',sCS.CSA.Car1(3),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(3) dOggParam(1) dOggParam(2)],...
                              'Visible',sCS.CSA.Car1(2),'Enable','on');
        oggIG.Car1 = [datoCar1 testoCar1];

        testoCar2 = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',sCS.CSA.Car2(1),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(4)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible',sCS.CSA.Car2(2),'Enable','on');
         datoCar2 = uicontrol('Style','edit',...
                              'String',sCS.CSA.Car2(3),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(4) dOggParam(1) dOggParam(2)],...
                              'Visible',sCS.CSA.Car2(2),'Enable','on');
        oggIG.Car2 = [datoCar2 testoCar2];

        testoCar3 = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',sCS.CSA.Car3(1),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(5)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible',sCS.CSA.Car3(2),'Enable','on');
         datoCar3 = uicontrol('Style','edit',...
                              'String',sCS.CSA.Car3(3),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(5) dOggParam(1) dOggParam(2)],...
                              'Visible',sCS.CSA.Car3(2),'Enable','on');
        oggIG.Car3 = [datoCar3 testoCar3];


            % % % Riga % % %
          TestoDt = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Dt,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(7)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
           datoDt = uicontrol('Style','edit',...
                              'String',sCS.CSA.DeltaT,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-3*(dOggParam(1)/2) rigRif(7) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.DeltaT = [datoDt TestoDt];

    testoNaturaP = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.NaturaP,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(7)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
     datoNaturaP = uicontrol('Style','popupmenu','Value',sCS.CSA.NaturaP,...
                              'FontSize',dimTOgg,'String',["Parabolica (∂/∂t)" "Iperbolica (∂²/∂t²)"],'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg+dOggParam(1)/2 rigRif(7) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','off');
        oggIG.NaturaP = [datoNaturaP testoNaturaP];

  testoIntervallo = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Intervallo,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(7)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
   datoIntervallo = uicontrol('Style','edit',...
                              'String',sCS.CSA.Intervallo,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(7) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.Intervallo = [datoIntervallo testoIntervallo];



        staccoOgg = 0.09;
        dOggTSpunta = [.08 .04]; % «[width height]»
            % % % Riga % % %
         datoStab = uicontrol('Style','checkbox','Value',sCS.CSA.Stabilita,...
                              'FontSize',dimTOgg,'String',testoOgg.Stab,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-staccoOgg-dOggTSpunta(1) rigRif(6)+staccoTesto/2 dOggTSpunta(1) dOggTSpunta(2)],...
                              'CallBack',@FDStabilita, ...      % colRif(1)-staccoOgg-dOggTSpunta(1)/2
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.Stabilita = datoStab;

        largdTTipoSim = 0.055; largDTipoSim = 0.07; distanzaTD = 0.0325;
     testoTipoSim = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','Sim. in','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-largdTTipoSim/2-distanzaTD rigRif(6)+staccoTesto/2 largdTTipoSim dOggParam(2)],...
                              'Visible','on','Enable','on');    % colRif(1)-staccoOgg-3*(dOggParam(1)/2) «[left bottom width height]»
     datoTipoSim = uicontrol('Style','popupmenu','Value',sCS.CSA.TipoSim,...
                              'String',["spazio" "tempo"],'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-largDTipoSim/2+distanzaTD rigRif(6)+staccoTesto/2 largDTipoSim dOggParam(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.TipoSim = [datoTipoSim testoTipoSim];


        datoMassa = uicontrol('Style','checkbox','Value',sCS.CSA.Ammassamento,...
                              'FontSize',dimTOgg,'String',testoOgg.Ammassamento,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)+staccoOgg rigRif(6)+staccoTesto/2 dOggTSpunta(1) dOggTSpunta(2)],...
                              'CallBack',@FDAmmassamento, ...   % colRif(1)+staccoOgg-dOggTSpunta(1)/2
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.Ammassamento = datoMassa;        


            % % % Terza colonna % % %
     testoAreaMax = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Area_Max,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(2)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
      datoAreaMax = uicontrol('Style','edit',...
                              'String',sCS.CSA.Area_Max,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(2) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.AreaMax = [datoAreaMax testoAreaMax];

   testoAngoloMin = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Theta_min,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(3)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
    datoAngoloMin = uicontrol('Style','edit',...
                              'String',sCS.CSA.Theta_min,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(3) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.AngoloMin = [datoAngoloMin testoAngoloMin];

        testoNSim = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.NSimulazioni,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(4)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
         datoNSim = uicontrol('Style','edit',...
                              'String',sCS.CSA.NSimulazioni,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(4) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.NSimulazioni = [datoNSim testoNSim];

         testoBase = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.Base,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(5)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
         datoBase = uicontrol('Style','edit',...
                              'String',sCS.CSA.Base,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggParam(1)/2 rigRif(5) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.Base = [datoBase testoBase];

                    uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.SelParam,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-(dOggFDSel(1)+0.01)/2 rigRif(1)+staccoTesto dOggFDSel(1)+0.01 dOggFDSel(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
                    uicontrol('Style','popupmenu',...
                              'FontSize',dimTOgg, ...
                              'String',sParam.NParam(2:end),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(1)-dOggFDSel(1)/2 rigRif(1) dOggFDSel(1) dOggFDSel(2)],...
                              'CallBack',@FDParametri, ...
                              'Visible','on','Enable','on');    % «[left bottom width height]»


            % % % Sesta colonna % % %
        testorifu = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(2)+staccoTesto dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
         datorifu = uicontrol('Style','edit',...
                              'String',sCS.CSA.Rifu,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(2) dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
        oggIG.Rifu = [datorifu testorifu];

      testodrifux = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(3)+staccoTesto dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
       datodrifux = uicontrol('Style','edit',...
                              'String',sCS.CSA.dRifu(1),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(3) dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
        oggIG.dRifux = [datodrifux testodrifux];

      testodrifuy = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(4)+staccoTesto dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
        datodrifuy = uicontrol('Style','edit',...
                              'String',sCS.CSA.dRifu(2),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(4) dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
        oggIG.dRifuy = [datodrifuy testodrifuy];

        testoriff = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(5)+staccoTesto dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
         datoriff = uicontrol('Style','edit',...
                              'String',sCS.CSA.Riff,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(5) dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
        oggIG.Riff = [datoriff testoriff];


        dOggTSpunta = [.1 .04]; % «[width height]»
        staccoOgg = 0.1;
            % % % Riga % % %
       datoEvoluz = uicontrol('Style','checkbox','Value',sCS.CSA.Evoluzione,...
                              'FontSize',dimTOgg,'String',testoOgg.Evoluzione,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-staccoOgg-dOggTSpunta(1) rigRif(6)+staccoTesto/2 dOggTSpunta(1) dOggTSpunta(2)],...
                              'CallBack',@FDEvoluzione, ...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
        oggIG.Evoluzione = datoEvoluz;

     testoIstante = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','Istante di val.','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggParam(1)/2 rigRif(6)+staccoTesto dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
      datoIstante = uicontrol('Style','edit',...
                              'String',sCS.CSA.Istante,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggParam(1)/2 rigRif(6) dOggParam(1) dOggParam(2)],...
                              'Visible','on','Enable','on');
        oggIG.Istante = [datoIstante testoIstante];

          datoVer = uicontrol('Style','checkbox','Value',sCS.CSA.Verifica,...
                              'FontSize',dimTOgg,'String',testoOgg.Verifica,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)+staccoOgg rigRif(6)+staccoTesto/2 dOggTSpunta(1) dOggTSpunta(2)],...
                              'CallBack',@FDVerifica, ...                        % -dOggTSpunta(1)*0.25
                              'Visible','on','Enable','on');
        oggIG.Verifica = datoVer;


        testodrifut = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String','∂u/∂t','HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(7)+staccoTesto dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
         datodrifut = uicontrol('Style','edit',...
                              'String',sCS.CSA.dRifu(3),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFunz(1)/2 rigRif(7) dOggFunz(1) dOggFunz(2)],...
                              'Visible','on','Enable','on');
        oggIG.dRifut = [datodrifut testodrifut];

                    uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.SelFunz,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFDSel(1)/2 rigRif(1)+staccoTesto dOggFDSel(1) dOggFDSel(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
                    uicontrol('Style','popupmenu',...
                              'FontSize',dimTOgg, ...
                              'String',sFunz.NFunz(:,2:end),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(3)-dOggFDSel(1)/2 rigRif(1) dOggFDSel(1) dOggFDSel(2)],...
                              'CallBack',@FDFunzioni, ...
                              'Visible','on','Enable','on');    % «[left bottom width height]»


            % % % Quinta colonna % % %
         testoDom = uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.SelDom,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(2)-dOggFDSel(1)/2 rigRif(1)+staccoTesto dOggFDSel(1) dOggFDSel(2)],...
                              'Visible','on','Enable','on');
          datoDom = uicontrol('Style','popupmenu','Value',sDom.(sCS.CSA.Nome_Dominio).Valore_Dominio,...
                              'FontSize',dimTOgg,'String',sDom.NDom(2:end),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(2)-dOggFDSel(1)/2 rigRif(1) dOggFDSel(1) dOggFDSel(2)],...
                              'CallBack',@FDDomini,...
                              'Visible','on','Enable','on');
        oggIG.Dominio = [datoDom testoDom];


            % Chiamata delle funzione di richiamo per impostare ulteriormente l'interfaccia d'avvio
        FDVerifica(datoVer,[])
        FDEvoluzione(datoEvoluz,[])
        FDStabilita(datoStab,[])
        FDAmmassamento(datoMassa,[])


%%

        % % % Tasti d'interazione: «avvio» e «riavvio» % % %

        
            % Parametri
        db   = [0.2 0.1]; % Dimensione bottoni
        pvB = 0.05;
        dimT = 12;        % Dimensione testo


            % Riavvio
        RiAvvioBott = uicontrol('Style','pushbutton', ...
                                'FontSize',dimT,'String','Riavvia', ...
                                'Units', 'Normalized','Position',[colRif(2)-db(1)/2 pvB db(1) db(2)], ...
                                'CallBack',@RiAvvio, ...
                                'Enable','off');

            % Avvio
                      uicontrol('Style','pushbutton',...
                                'FontSize',dimT,'String','Avvia', ...
                                'Units', 'Normalized','Position',[colRif(1)-db(1)/2 pvB db(1) db(2)],...
                                'CallBack',{@Avvio, RiAvvioBott},...
                                'Enable','on');

            % Terminazione
                      uicontrol('Style','pushbutton',...
                                'FontSize',dimT,'String','Esci', ...
                                'Units', 'Normalized','Position',[colRif(3)-db(1)/2 pvB db(1) db(2)],...
                                'CallBack',@Esci,...
                                'Enable','on');

            % % % Quinta colonna % % %
        staccoSelCasi = 0.125;
                    uicontrol('Style','text',...
                              'FontSize',dimTOgg,'String',testoOgg.CS,'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(2)-dOggFDSel(1)/2 pvB+staccoTesto+staccoSelCasi dOggFDSel(1) dOggFDSel(2)],...
                              'Visible','on','Enable','on');    % «[left bottom width height]»
                    uicontrol('Style','popupmenu',...
                              'FontSize',dimTOgg, ...
                              'String',sCS.NCS(2:end),'HorizontalAlignment','center',...
                              'Units', 'Normalized','Position',[colRif(2)-dOggFDSel(1)/2 pvB+staccoSelCasi dOggFDSel(1) dOggFDSel(2)],...
                              'CallBack',@FDCasiStudio, ...
                              'Visible','on','Enable','on');    % «[left bottom width height]»


        % Ciclo «while» per fermare il programma finché il bottone «AvvioBott» non è stato premuto
        while(avvia == 0 && termina == 0)
            pause(0.01) % Pausa di 0.01 secondi: necessaria per dare respiro a Matlab, altrimenti neanche genera la figura precedente
        end

%%
        
            %%% Traduzione delle stringhe inserite in tutte i dati necessari alla simulazione %%%


            datiS = []; % Dati sulla simulazione
        datiS.ordP      = datoOrd.Value;
        datiS.verifica  = datoVer.Value;

        datiS.stabilita = datoStab.Value;
        datiS.massa     = datoMassa.Value;

        datiS.angoloMin = str2double(datoAngoloMin.String);
        datiS.areaMax   = str2double(datoAreaMax.String);
        datiS.base      = str2double(datoBase.String);

        datiS.nSim      = str2double(datoNSim.String);
        datiS.tipoSim   = datoTipoSim.Value;

        
            datiE = []; % Dati sull'equazione
        datiE.ni    = str2num(datoNi.String,Evaluation="restricted");
        datiE.beta  = str2num(datoBeta.String,Evaluation="restricted")';
        datiE.sigma = str2num(datoSigma.String,Evaluation="restricted");
        

            datiD = []; % Dati sul dominio
        switch datoDom.Value
            case 1
                datiD.tipoDom = "quadrato";
                datiD.infoDom = [str2double(datoCar1.String) str2num(datoCar2.String,Evaluation="restricted")];
                datiD.lati    = str2num(datoLati.String,Evaluation="restricted");
            case 2
                datiD.tipoDom = "elle";
                datiD.infoDom = [str2double(datoCar1.String) str2num(datoCar2.String,Evaluation="restricted")];
                datiD.lati    = str2num(datoLati.String,Evaluation="restricted");
            case 3
                datiD.tipoDom = "triangolo";
                datiD.infoDom = [str2double(datoCar1.String) ...
                                 str2double(datoCar2.String) ...
                                 str2num(datoCar3.String,Evaluation="restricted")];
                datiD.lati    = str2num(datoLati.String,Evaluation="restricted");
            case 4
                datiD.tipoDom = "cerchio";
                datiD.infoDom = [str2double(datoCar1.String) ...
                                 str2num(datoCar2.String,Evaluation="restricted") ...
                                 str2num(datoCar3.String,Evaluation="restricted")];
                
                datiD.nDir = floor(str2double(datoLati.String)*str2double(datoCar1.String));
                datiD.lati = [ones(1,nDir)  2*ones(1,infoDom(1)-nDir)]; % Dirichlet<->1 e Neumann<->2
        end

            datiF = []; % Dati sulle funzioni
        datiF.rifu   = str2func(append('@(x,y,t)',datorifu.String));
        datiF.drifux = str2func(append('@(x,y,t)',datodrifux.String));
        datiF.drifuy = str2func(append('@(x,y,t)',datodrifuy.String));
        datiF.riff   = str2func(append('@(x,y,t)',datoriff.String));
            % https://it.mathworks.com/help/matlab/ref/str2func.html

            datiT = []; % Dati sull'evoluzione
        datiT.evol    = datoEvoluz.Value;

        datiT.dt      = str2num(datoDt.String,Evaluation="restricted");
        datiT.T       = str2num(datoIntervallo.String,Evaluation="restricted");
        datiT.natura  = datoNaturaP.Value;
        datiT.I       = str2num(datoIstante.String,Evaluation="restricted");

        datiT.drifut   = str2func(append('@(x,y,t)',datodrifut.String));


end


% Le funzioni di richiamo («callback») richiedono due ingressi: «source» e «event»; la prima contiene le informazioni
% su ciò che è stato premuto d'altra parte mentre la seconda è inutile per cui è sufficiente ignorarli mediante «~»

%%

% Per chiamare una funzione di richiamo dentro un'altra funzione di richiamo si veda la preziosa risposta qui
% (https://stackoverflow.com/questions/1819232/how-do-i-execute-a-callback-function-from-another-function-file-in-matlab)

function FDDomini(srg,~) % F[inestra di ]D[ialogo dei ]Domini

    global oggIG sDom

    campiDom = string(fieldnames(sDom));

        % Si rendono invisibili tutte le precedenti finestre [per tenere conto di precedenti selezioni]
    oggIG.Car1(2).Visible = "off"; oggIG.Car1(1).Visible = "off";
    oggIG.Car2(2).Visible = "off"; oggIG.Car2(1).Visible = "off";
    oggIG.Car3(2).Visible = "off"; oggIG.Car3(1).Visible = "off";
        % (1) sta per il dato e (2) per testo e cosí vale per tutte le chiamate successive


        % Si cambia il testo e s'impone la percentuale di lati di Dirichlet uguale a 1 (tutti Dirichlet)
    oggIG.Lati(2).String = sDom.(campiDom(srg.Value+1)).Info_Contorno(1);
    oggIG.Lati(1).String = sDom.(campiDom(srg.Value+1)).Info_Contorno(2); % Dirichlet<->1 e Neumann<->2
    
        % Si rendono visibili le prime tre finestre secondo le caratteristiche del caso di studio
    oggIG.Car1(2).Visible = sDom.(campiDom(srg.Value+1)).Car1(2); oggIG.Car1(1).Visible = sDom.(campiDom(srg.Value+1)).Car1(2);
    oggIG.Car2(2).Visible = sDom.(campiDom(srg.Value+1)).Car2(2); oggIG.Car2(1).Visible = sDom.(campiDom(srg.Value+1)).Car2(2);
    oggIG.Car3(2).Visible = sDom.(campiDom(srg.Value+1)).Car3(2); oggIG.Car3(1).Visible = sDom.(campiDom(srg.Value+1)).Car3(2);

        % Si rinominano le stringhe delle caratteristiche
    oggIG.Car1(2).String = sDom.(campiDom(srg.Value+1)).Car1(1); oggIG.Car1(1).String = sDom.(campiDom(srg.Value+1)).Car1(3);
    oggIG.Car2(2).String = sDom.(campiDom(srg.Value+1)).Car2(1); oggIG.Car2(1).String = sDom.(campiDom(srg.Value+1)).Car2(3);
    oggIG.Car3(2).String = sDom.(campiDom(srg.Value+1)).Car3(1); oggIG.Car3(1).String = sDom.(campiDom(srg.Value+1)).Car3(3);

end

function FDParametri(srg,~) % F[inestra di ]D[ialogo dei ]Parametri

    global oggIG sParam

    campiParam = string(fieldnames(sParam));

    oggIG.Ni(1).String    = sParam.(campiParam(srg.Value+1)).Ni;
    oggIG.Beta(1).String  = sParam.(campiParam(srg.Value+1)).Beta;
    oggIG.Sigma(1).String = sParam.(campiParam(srg.Value+1)).Sigma;

end

function FDVerifica(srg,~) % F[inestra di ]D[ialogo della ]Verifica

    global oggIG

    if(srg.Value == 1) % Se v'è la spunta
            % Le tre caselle corrispondono alle vere funzioni da cui discendono gD e gN = [gNx gNy]
        oggIG.Rifu(2).String   = "Funzione u";
        oggIG.dRifux(2).String = "Funzione ∂u/∂x";
        oggIG.dRifuy(2).String = "Funzione ∂u/∂y";
        oggIG.Riff(2).String   = "Funzione ∆u";
        oggIG.dRifut(2).String = "Funzione ∂u/∂t";
    else               % Altrimenti
            % Le tre caselle corrispondono direttamente alle funzioni gD e gN = [gNx gNy]
        oggIG.Rifu(2).String   = "Funzione gD";
        oggIG.dRifux(2).String = "Funzione gNx";
        oggIG.dRifuy(2).String = "Funzione gNy";
        oggIG.Riff(2).String   = "Funzione f";
        oggIG.dRifut(2).String = "Funzione u0";
    end

end

function FDFunzioni(srg,~) % F[inestra di ]D[ialogo delle ]Funzioni

    global oggIG sFunz

    campiFunz = string(fieldnames(sFunz));

        % Secondo «srg.Value» si seleziona se svolgere la verifica, e quindi le tre caselle corrispondono
        % alle vere funzioni da cui discendono gD e gN = [gNx gNy], oppure no, e quindi le tre caselle
        % corrispondono direttamente alle funzioni gD e gN = [gNx gNy]; ciò viene fatto chiamando FDVerifica
    oggIG.Verifica.Value = sFunz.(campiFunz(srg.Value+1)).Verifica;
    FDVerifica(oggIG.Verifica)

        % S'impostano le funzioni relativo al caso identificato da «srg.Value»
    oggIG.Rifu(1).String   = sFunz.(campiFunz(srg.Value+1)).Rifu;
    oggIG.dRifux(1).String = sFunz.(campiFunz(srg.Value+1)).dRifu(1);
    oggIG.dRifuy(1).String = sFunz.(campiFunz(srg.Value+1)).dRifu(2);
    oggIG.Riff(1).String   = sFunz.(campiFunz(srg.Value+1)).Riff;

end

function FDCasiStudio(srg,~) % F[inestra di ]D[ialogo dei ]Casi[ di ]Studio

    global oggIG sCS sDom

    campiDom = string(fieldnames(sCS));

        % 1° Pezzo
    oggIG.Verifica(1).Value = sCS.(campiDom(srg.Value+2)).Verifica;
    FDVerifica(oggIG.Verifica(1))

    oggIG.OrdineP(1).Value      = sCS.(campiDom(srg.Value+2)).OrdineP;
    oggIG.Stabilita(1).Value    = sCS.(campiDom(srg.Value+2)).Stabilita;
    oggIG.Ammassamento(1).Value = sCS.(campiDom(srg.Value+2)).Ammassamento;
    FDStabilita(oggIG.Stabilita(1),[])
    FDAmmassamento(oggIG.Ammassamento(1),[])

    oggIG.Ni(1).String    = sCS.(campiDom(srg.Value+2)).Ni;
    oggIG.Beta(1).String  = sCS.(campiDom(srg.Value+2)).Beta;
    oggIG.Sigma(1).String = sCS.(campiDom(srg.Value+2)).Sigma;

        % 2° Pezzo
    oggIG.AreaMax(1).String      = sCS.(campiDom(srg.Value+2)).Area_Max;
    oggIG.AngoloMin(1).String    = sCS.(campiDom(srg.Value+2)).Theta_min;
    oggIG.NSimulazioni(1).String = sCS.(campiDom(srg.Value+2)).NSimulazioni;
    oggIG.Base(1).String         = sCS.(campiDom(srg.Value+2)).Base;

        % 3°, 4° e 5° Pezzo
    oggIG.Car1(2).Visible = "off"; oggIG.Car1(1).Visible = "off"; 
    oggIG.Car2(2).Visible = "off"; oggIG.Car2(1).Visible = "off";
    oggIG.Car3(2).Visible = "off"; oggIG.Car3(1).Visible = "off";

        % Si cambia la geometria, il testoCaratt e s'impone la percentuale di lati di Dirichlet uguale a 1 (tutti Dirichlet)
    oggIG.Dominio(1).Value = sDom.(sCS.(campiDom(srg.Value+1)).Nome_Dominio).Valore_Dominio;
    oggIG.Lati(2).String   = sCS.(campiDom(srg.Value+2)).Info_Contorno(1); 
    oggIG.Lati(1).String   = sCS.(campiDom(srg.Value+2)).Info_Contorno(2); % Dirichlet<->1 e Neumann<->2
    
        % Si rendono visibili le prime tre finestre
    oggIG.Car1(2).Visible = sCS.(campiDom(srg.Value+2)).Car1(2); oggIG.Car1(1).Visible = sCS.(campiDom(srg.Value+2)).Car1(2);
    oggIG.Car2(2).Visible = sCS.(campiDom(srg.Value+2)).Car2(2); oggIG.Car2(1).Visible = sCS.(campiDom(srg.Value+2)).Car2(2);
    oggIG.Car3(2).Visible = sCS.(campiDom(srg.Value+2)).Car3(2); oggIG.Car3(1).Visible = sCS.(campiDom(srg.Value+2)).Car3(2);

        % Si rinominano le stringhe delle caratteristiche
    oggIG.Car1(2).String = sCS.(campiDom(srg.Value+2)).Car1(1); oggIG.Car1(1).String = sCS.(campiDom(srg.Value+2)).Car1(3);
    oggIG.Car2(2).String = sCS.(campiDom(srg.Value+2)).Car2(1); oggIG.Car2(1).String = sCS.(campiDom(srg.Value+2)).Car2(3);
    oggIG.Car3(2).String = sCS.(campiDom(srg.Value+2)).Car3(1); oggIG.Car3(1).String = sCS.(campiDom(srg.Value+2)).Car3(3);

    oggIG.Rifu(1).String   = sCS.(campiDom(srg.Value+2)).Rifu;
    oggIG.dRifux(1).String = sCS.(campiDom(srg.Value+2)).dRifu(1);
    oggIG.dRifuy(1).String = sCS.(campiDom(srg.Value+2)).dRifu(2);
    oggIG.Riff(1).String   = sCS.(campiDom(srg.Value+2)).Riff;

        % Tempo
    oggIG.DeltaT(1).String      = sCS.(campiDom(srg.Value+2)).DeltaT;
    oggIG.Intervallo(1).String  = sCS.(campiDom(srg.Value+2)).Intervallo;
    oggIG.NaturaP(1).Value  = sCS.(campiDom(srg.Value+2)).NaturaP;
    oggIG.dRifut(1).String   = sCS.(campiDom(srg.Value+2)).dRifu(3);
    oggIG.Istante(1).String = sCS.(campiDom(srg.Value+2)).Istante;

    oggIG.TipoSim(1).Value = sCS.(campiDom(srg.Value+2)).TipoSim;

    oggIG.Evoluzione.Value = sCS.(campiDom(srg.Value+2)).Evoluzione;
    FDEvoluzione(oggIG.Evoluzione)

end

function Avvio(srg,~,RiAvvioBott)
        % Variabile globale per tracciare se è stato premuto il bottone; infatti Matlab, per qualche motivazione,
        % cambia «AvvioBott.Value» da 0 a 1 solo finché si è all'interno della funzione di richiamo: una volta fuori il
        % suo valore è riportato a zero, come se non fosse stato premuto
    global avvia
    avvia = 1;

        % Disabilità il bottone d'avvio e abilità quello di riavvio
    srg.Enable         = 'off';
    RiAvvioBott.Enable = 'on';
end

function RiAvvio(~,~)
    run Es_P1P2_EqCalore.m;
end

function Esci(~,~)
    global termina avvia

    if(avvia == 0) % Se non si è avviato il codice
        termina = 1;
    else           % Se non si è avviato il codice
        close all
    end
end

function FDStabilita(srg,~)
    global oggIG
    if(srg.Value == 1)
        oggIG.Ammassamento.Value = 0;
        oggIG.Ammassamento.Enable = "off";
    else
        oggIG.Ammassamento.Enable = "on";
    end        
end

function FDAmmassamento(srg,~)
    global oggIG
    if(srg.Value == 1)
        oggIG.Stabilita.Value = 0;
        oggIG.Stabilita.Enable = "off";
    else
        oggIG.Stabilita.Enable = "on";
    end        
end

function FDEvoluzione(srg,~)
    
    global oggIG

    if(srg.Value == 0)

        oggIG.DeltaT(1).Enable      = "off";
        oggIG.Intervallo(1).Enable  = "off";
        oggIG.Intervallo(1).String  = "0";
        oggIG.Istante(2).String     = "Istante di val.";
        oggIG.TipoSim(1).Value         = 1;
        oggIG.TipoSim(1).Enable        = "off";
    else
        oggIG.DeltaT(1).Enable      = "on";
        oggIG.Intervallo(1).Enable  = "on";
        oggIG.Istante(2).String     = "Istante iniziale";
        oggIG.TipoSim(1).Enable        = "on";
    end

end