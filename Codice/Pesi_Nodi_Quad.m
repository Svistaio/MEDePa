function [xhat, yhat, omega] = Pesi_Nodi_Quad()

% Questa funzione approssima esattamente polinomi di grado cinque che va abbastanza bene fino agli elementi finiti P2; infatti nell'equazione di riferimento il prodotto tra funzioni sarà al massimmo, appunto, di grado cinque se si hanno i P2

sqrt15 = sqrt(15.0);

xhat = [ 1/3
	    (6 + sqrt15) / 21
	    (9 - 2*sqrt15) / 21
	    (6 + sqrt15) / 21
	    (6 - sqrt15) / 21
	    (9 + 2*sqrt15) / 21
	    (6 - sqrt15) / 21];

yhat = [1/3
        (6 + sqrt15) / 21
	    (6 + sqrt15) / 21
	    (9 - 2*sqrt15) / 21
	    (6 - sqrt15) / 21
	    (6 - sqrt15) / 21
	    (9 + 2*sqrt15) / 21];

% sum(omega)=1/2
% int_T 1 d Omega = |T|
omega = [ 9/80
          (155 + sqrt15) / 2400
          (155 + sqrt15) / 2400
          (155 + sqrt15) / 2400
          (155 - sqrt15) / 2400
          (155 - sqrt15) / 2400
          (155 - sqrt15) / 2400];

return
