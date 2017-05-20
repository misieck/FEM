    Q = 100;
    A = 10;
    k = 5;
    Nelem=160;
    Ndof = Nelem +1;
    L  = 6/Nelem;
    T0 = 100;
    q_4_boundary = 15;
    bc = [1 T0]
    
    Edof = zeros(Nelem, 3);
    
    for i=1:Nelem 
        Edof(i,:) = [i, i, i + 1];
    end
       
   EssentialNodeNo = 1;
    
   %materialegenskapsmatris
   K = zeros(Ndof);
   
   %lastvektorn
   fl = zeros(Ndof, 1);
   
   %elementmaterialegenskap
   ep = k*A/L %simple
   Ke = spring1e(ep) %matrix
   
   %elementlastvektorn
   fe = [100*L/2; 100*L/2];
   
for n = 1:Nelem
   [K, fl]=assem(Edof(n,:),K,Ke, fl, fe);
end


%cut out essential node from equation
fb_clipped = [ zeros(Ndof-2, 1); -A*q_4_boundary];

K_clipped = K;
K_clipped (EssentialNodeNo,:) = [];
K_clipped (:, EssentialNodeNo) = [];
fl_clipped = fl;
fl_clipped(EssentialNodeNo) = [];
f_clipped = fl_clipped + fb_clipped;




%Denna borde kallas A
%[T_clipped, Q]=solveq(K_clipped,f_clipped)
%T = [T0;T_clipped]
[T, Q]=solveq(K,[100*L/2;f_clipped], bc);

Q(1)
%q_1_boundary = spring1s( ep, [T(2), T(1)])

T(1:4)