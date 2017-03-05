# MULPHYS
cd /dat/mulphys/run; 
./mulphys -f TAM.xml # t T b B g G
./mulphys -f ParticleFlow.xml # gui.cfg: Vector.length=2 l: v 5 0 G 
./mulphys -f Membrane.xml # f GG B gui.cfg: Vector.length=8 l: v 10 0 B b
./mulphys -f Lattice.xml # n, e
cd elastic; eog -f *.gif
cd dove; eog -f *.png *.gif &
cd collapse; eog -f *.png *.gif &
cd biomed; eog -f *.png *.gif &
cd nano; eog -f *.png *.gif &
# RASAT
cd /dat/www/html/mulphys/rasat/demo; appletviewer city.html
cd /dat/www/html/mulphys/mesh/mc/demo/run; appletviewer square.html
cd /dat/www/html/mulphys/mesh/mc/demo/run; appletviewer cube.html
# REMODY
cd /dat/remody
