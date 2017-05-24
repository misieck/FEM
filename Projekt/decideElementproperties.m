

function [E,nu,k,ro,c,alpha] = decideElementproperties(elements, nelem)
     constants;
    
     E = zeros(nelem, 1);
     nu = zeros(nelem, 1);
     k = zeros(nelem, 1);
     ro = zeros(nelem, 1);
     c = zeros(nelem, 1);
     alpha = zeros(nelem, 1);
    
    for i = 1:nelem
     switch elements(4, i)
        case 1
            
            E(i) = E_pcb;
            nu(i) = nu_pcb;
            k(i) = k_pcb;
            ro(i) = ro_pcb;
            c(i) = c_pcb;
            alpha(i) = alpha_pcb;
            
        case 2
            
            E(i) = E_smd;
            nu(i) = nu_smd;
            k(i) = k_smd;
            ro(i) = ro_smd;
            c(i) = c_smd;
            alpha(i) = alpha_smd;
         
        case 3
            
            E(i) = E_sol;
            nu(i) = nu_sol;
            k(i) = k_sol;
            ro(i) = ro_sol;
            c(i) = c_sol;
            alpha(i) = alpha_sol;
     end
    end

    
