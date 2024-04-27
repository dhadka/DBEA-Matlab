function ind = realmutate(a,rate,domains)
%REALMUTATE Summary of this function goes here
%   Detailed explanation goes here

   % double rnd, delta1, delta2, mut_pow, deltaq;
   % double y, yl, yu, val, xy;
   % double eta_m = id_mu;
   nx=size(a,2);
   eta_m=20;
   for j = 1:nx
        if (rand() <= rate) 
            y = a(j);
            
            yl = domains(j,1);
            yu = domains(j,2);
            delta1 = (y - yl) / (yu - yl);
            delta2 = (yu - y) / (yu - yl);

            rnd = rand();
            mut_pow = 1.0 / (eta_m + 1.0);
            if (rnd <= 0.5) 
                xy = 1.0 - delta1;
                val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (xy^(eta_m + 1.0));
                deltaq = (val^mut_pow) - 1.0;
            else 
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (xy^ (eta_m + 1.0));
                deltaq = 1.0 - (val^mut_pow);
            end
            
            y = y + deltaq * (yu - yl);
            if (y < yl)
                y = yl;
            end
            if (y > yu)
                y = yu;
            end
            a(j) = y;        
        end
   end
  ind = a;
end