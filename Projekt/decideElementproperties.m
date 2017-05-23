

function [E,nu,k,ro,c,alpha,D] = decideElementproperties(elements, IndexForElement)
     constants;
        
     switch elements(4, IndexForElement)
        case 1
            
            E = E_pcb;
            nu = nu_pcb;
            k = k_pcb;
            ro = ro_pcb;
            c = c_pcb;
            alpha = alpha_pcb;
            
            D = D_pcb;
            
        case 2
            
            E = E_smd;
            nu = nu_smd;
            k = k_smd;
            ro = ro_smd;
            c = c_smd;
            alpha = alpha_smd;
            
            D = D_smd;
            
        case 3
            
            E = E_sol;
            nu = nu_sol;
            k = k_sol;
            ro = ro_sol;
            c = c_sol;
            alpha = alpha_sol;
            
            D = D_sol;
     end

    
