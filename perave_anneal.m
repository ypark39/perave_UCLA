f=@(y) perave_opti(y(1),y(2),-20,.66);
simulannealbnd(f,[16e-6,1.8],[5e-6,0],[20e-6,2*pi]);