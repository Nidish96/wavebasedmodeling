function [] = PLOTSOLN(U,BM1, BM2, IN1, IN2, L1, L2, Nemono, Nein, sc, ccm1, ccm2, alph)
  plot(BM1.X, BM1.Y, 'ko-'); hold on 
  plot(IN1.X, IN1.Y, 'ko-'); hold on
  plot(BM2.X, BM2.Y, 'ko-'); hold on 
  plot(IN2.X, IN2.Y, 'ko-'); hold on 
  grid on 

  PLANARBMDEPICT(L1(1:(Nemono+1)*3, :)*U*sc, BM1, ccm1, 0.6);
  PLANARBMDEPICT(L1((Nemono+1)*3+(1:(Nein+1)*3), :)*U*sc, IN1, ccm1, 0.6);

  PLANARBMDEPICT(L2((Nein+1)*3+(1:(Nemono+1)*3), :)*U*sc, BM2, ccm2, 0.6);
  PLANARBMDEPICT(L2(1:(Nein+1)*3, :)*U*sc, IN2, ccm2, 0.6);
end
