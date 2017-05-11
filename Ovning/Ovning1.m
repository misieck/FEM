    Q = 100;
    A = 10;
    k = 5;

    Edof = [1 1 2;
        2 2 3;
        3 3 4];
   
   K = zeros(4);
   f = zeros(4, 1);
   
   ep = k*A;
   
   Ke = spring1e(ep);
   fe = [100; 100];
   
for n = 1:3
   [f, K]=assem(Edof(n,:),K,Ke, f, fe);
end

fb = 